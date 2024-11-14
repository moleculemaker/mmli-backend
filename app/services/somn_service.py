import io
import os
import json
import uuid

from fastapi import HTTPException
import pandas as pd
from sqlmodel.ext.asyncio.session import AsyncSession

from config import get_logger
from models.sqlmodel.models import Job

from services.minio_service import MinIOService

from typing import List, Literal
from rdkit import Chem
from openbabel import pybel as pb
import traceback

log = get_logger(__name__)
from openbabel import openbabel as ob

SOMN_ERROR_TYPES = Literal[
    'invalid_input', 
    '3d_gen', 
    'no_reactive_nitrogens', 
    'no_br_or_cl_in_el', 
    'no_nitrogens_in_nuc',
    'br_in_nuc',
    'cl_in_nuc',
    'no_reaction_site'
]

class SomnException(BaseException):
    type: SOMN_ERROR_TYPES
    message: str
    
    def __init__(self, type: str, message: str):
        self.type = type
        self.message = message

class SomnService:
    somn_frontend_baseURL = os.environ.get("SOMN_FRONTEND_URL")
    
    @staticmethod
    def has_chiral(mol: pb.Molecule):
        '''
        Check if the molecule has any chiral centers (so we can prompt warning/info message) on front-end

        Parameters:
            mol (class 'openbabel.pybel.Molecule') - pybel.molecule object read from the SMILES String
        
        Returns:
            Bool - True: has stereochemistry, False: no stereochemistry 
        '''
        m = mol.OBMol
        for genericdata in m.GetAllData(ob.StereoData):
            stereodata = ob.toStereoBase(genericdata)
            stereotype = stereodata.GetType()
            if (stereotype):
                return True
        return False
    
    @staticmethod
    def get_num_heavy_atoms(mol: pb.Molecule):
        """
        Count the number of (non-Hydrogen) heavy atoms

        Params:
            mol (class 'openbabel.pybel.Molecule') - pybel.molecule object read from the SMILES String 

        Returns:
            cntHeavyAtoms (int) - number of heavy atoms in the molecule
        """

        # Remove Hs
        mol.removeh()

        # Count atoms
        myVec=mol.atoms

        # Get number of non-H atoms
        cntHeavyAtoms=len(myVec)

        # Add back Hs for further processing
        mol.addh()

        return cntHeavyAtoms
    
    @staticmethod
    def canonicalize_smiles(smiles: str):
        mol = pb.readstring('smi', smiles)
        return mol.write('smi').strip()
    
    @staticmethod
    def ob_test(
        user_input: str,
        input_type: Literal['smi', 'cml', 'cdxml']
    ):
        """
        Generate 3D coordinates for a molecule from its SMILES string using OpenBabel.
        
        Parameters:
            smiles (str): SMILES representation of the molecule
            
        Returns:
            tuple[str, bool, int]: A tuple containing:
                - The molecule in MOL2 format with 3D coordinates
                - Boolean indicating if the molecule has chiral centers
                - Number of heavy atoms
                
        Raises:
            SomnException: If SMILES string is invalid or 3D coordinate generation fails
        """
        try:
            obmol = pb.readstring(input_type, user_input)
        except Exception as e:
            raise SomnException(type="invalid_input", message=f"Invalid input [{input_type}]: {user_input}")

        obmol.addh()
        
        try:            
            # this step may fail, so we know SOMN cannot compute on the input
            obmol.make3D() 
        except Exception as e:
            log.error(f"Unable to generate 3D coordinates {traceback.print_exc()}")
            raise SomnException(type="3d_gen", message=f"Unable to generate 3D coordinates for {user_input}")
        
        return (
            obmol.write("mol2"), 
            SomnService.has_chiral(obmol), 
            SomnService.get_num_heavy_atoms(obmol)
        )
    
    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        job = await db.get(Job, job_id)
        job_info = json.loads(job.job_info)
        ret_val = []
        for info in job_info['info']:
            file_name = f"{info['nuc_name']}_{info['el_name']}_processed.csv"
            
            csv_content = service.get_file(bucket_name, f"/{job_id}/out/{file_name}")
            if csv_content is None:
                raise HTTPException(status_code=404, detail=f"File {file_name} not found")
            
            df = pd.read_csv(io.BytesIO(csv_content))
            file_data = []
            for _, row in df.iterrows():
                data = {}
                keys = row[0].split('_')
                average = row['average']
                stdev = row['stdev']
                data['nuc_name'] = keys[0]
                data['el_name'] = keys[1]
                data['catalyst'] = int(keys[2])
                data['solvent'] = int(keys[3])
                data['base'] = keys[4]
                data['yield'] = average
                data['stdev'] = stdev
                file_data.append(data)
            ret_val.append(file_data)
            
        return ret_val
    
    @staticmethod
    def check_user_input_substrates(
        user_input, 
        input_type: Literal['smi', 'cml', 'cdxml'], 
        role: Literal['el', 'nuc']
    ):
        """
        Verifies user input substrate and returns reaction site indices and chirality information.

        Parameters:
            user_input (str): SMILES string of the input molecule
            role (str): Role of the molecule - either 'el' for electrophile or 'nuc' for nucleophile

        Returns:
            tuple: A tuple containing:
                - list[int]: List of atom indices for reaction sites
                - a mol2 string of the molecule with 3D coordinates
                - bool: True if molecule has chiral centers, False otherwise

        Raises:
            SomnException: If no valid reaction sites are found or if invalid molecule type
            Exception: For other validation errors
        """

        def get_amine_ref_ns(mol, ref_atom_idxes: List) -> List:
            """
            Returns the reference atom index for the nitrogen most likely to be reactive, or multiple (if there are multiple)
            """
            reactive_nitrogens = []
            reactive_nitrogen_idxes = []
            for idx in ref_atom_idxes:
                neighbors = mol.GetAtomWithIdx(idx).GetNeighbors()
                n_hs = [1 for nbr in neighbors if nbr.GetSymbol() == "H"]
                reactive_nitrogens.append(sum(n_hs))
                reactive_nitrogen_idxes.append(idx)

            reactive_nitrogens = np.array(reactive_nitrogens)
            atm_indices = np.where((reactive_nitrogens == np.max(reactive_nitrogens)) & (reactive_nitrogens != 0))[0]
            
            ret_val = []
            for i in atm_indices:
                ret_val.append(reactive_nitrogen_idxes[i])

            if len(ret_val) >= 1:
                return ret_val
            
            raise SomnException(type="no_reactive_nitrogens", message="No reactive nitrogens found")
            
        def check_halides_aromatic(rdkmol,halides):
            rdkatoms = [atom for atom in rdkmol.GetAtoms() if atom.GetIdx() in halides]

            outIdxes = []
            for atom in rdkatoms:
                nbrs = atom.GetNeighbors()
                assert len(nbrs) == 1
                if nbrs[0].GetIsAromatic() is True: 
                    outIdxes.append(atom.GetIdx())

            return outIdxes

        def get_atoms_by_symbol(mol, symbol):
            num_atoms = mol.GetNumAtoms()
            retVals = []
            for idx in range(num_atoms):
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == symbol:
                    retVals.append(idx)
            return retVals

        # add ob_test here to prevent the user from 
        # submitting a molecule that cannot be processed by SOMN
        (mol2, has_chiral, num_heavy_atoms) = SomnService.ob_test(user_input, input_type)

        # generate rdkit mol using mol2 because rdkit
        # might have issues generating mol from smiles directly
        mol = Chem.MolFromMol2Block(mol2, sanitize=False, removeHs=False)

        bromides = get_atoms_by_symbol(mol, symbol="Br")
        chlorides = get_atoms_by_symbol(mol, symbol="Cl")
        nitrogens = get_atoms_by_symbol(mol, symbol="N")

        if role.startswith('nuc'):
            import numpy as np
            if len(bromides) != 0:
                aromatic_halides = check_halides_aromatic(mol,bromides)
                if any(idx in aromatic_halides for idx in bromides):
                    raise SomnException(type="br_in_nuc", message="Bromine found in nucleophile")

            if len(chlorides) != 0:
                aromatic_halides = check_halides_aromatic(mol,bromides)
                if any(idx in aromatic_halides for idx in chlorides):
                    raise SomnException(type="cl_in_nuc", message="Chlorine found in nucleophile")
            
            if len(nitrogens) == 0:
                raise SomnException(type="no_nitrogens_in_nuc", message="No nitrogens found in nucleophile")
            
            indices = get_amine_ref_ns(mol,nitrogens)
            if not len(indices):
                raise SomnException(type="no_reaction_site", message="No Reaction Site found")
            
            return (indices, mol2, has_chiral, num_heavy_atoms)

        elif role.startswith('el'):
            if len(bromides) + len(chlorides) == 0:
                raise SomnException(type="no_br_or_cl_in_el", message="No Br or Cl found in electrophile")
            
            try:
                chl_idxes = check_halides_aromatic(mol,chlorides)
                brm_idxes = check_halides_aromatic(mol,bromides)
                
                if not len(chl_idxes) and not len(brm_idxes):
                    raise SomnException(type="no_reaction_site", message="No Reaction Site found")

                return (chl_idxes + brm_idxes, mol2, has_chiral, num_heavy_atoms)
            
            except Exception as e:
                raise e

        raise SomnException(type="invalid_input", message=f"Invalid input role: {role}")
    
    @staticmethod
    def process_molecule_input(config: dict, mol_type: Literal['el', 'nuc']) -> str:
        """Process molecule input and return standardized format."""
        reaction_sites, mol_str, _, _ = SomnService.check_user_input_substrates(
            config[mol_type], 
            config[f'{mol_type}_input_type'], 
            mol_type
        )
        
        if len(reaction_sites) > 1 or config[f'{mol_type}_input_type'] in ['cdxml', 'cml']:
            return mol_str
        return SomnService.canonicalize_smiles(config[mol_type])
    
    @staticmethod
    def generate_name_mapping(configs: List[dict], mol_type: Literal['el', 'nuc']) -> dict:
        """Generate unique name mappings for molecules."""
        name_map = {}
        for config in configs:
            name_key = f'{mol_type}_name'
            if config[name_key] not in name_map:
                new_name = str(uuid.uuid4()).replace('-', '')
                name_map[config[name_key]] = new_name
            config[name_key] = name_map[config[name_key]]
        return (configs, name_map)
    
    @staticmethod
    def update_names_from_reference(
        configs: List[dict], 
        name_map: dict,
        mol_type: Literal['el', 'nuc']
    ) -> None:
        """Update names from reference JSON file."""
        ref = {
            'el': {
                "Brc1c(C(C)C)cc(C(C)C)cc1C(C)C": "1",
                "Brc1c(CC)cccc1CC": "10",
                "Brc1cc(SC)ccc1": "11",
                "Brc1c2c(n(C)nc2)ccc1": "12",
                "Brc1cc([N](=O)[O-])ccc1": "13",
                "O=c1n(C)cnc2c1cc(Br)cc2": "14",
                "Brc1ccccc1": "15",
                "Brc1cc(O[Si](C(C)C)(C(C)C)C(C)C)ccc1": "16",
                "Brc1cnc(C)nc1": "17",
                "Brc1cnc(C)cc1": "18",
                "Brc1nccs1": "19",
                "Brc1c(OC)cccc1OC": "2",
                "Brc1cscn1": "20",
                "Brc1cocc1": "21",
                "Brc1cc2ccccc2nc1": "22",
                "Brc1cc(OC)c(OC)c(OC)c1": "23",
                "Brc1cncnc1": "24",
                "Brc1cc2c(nc(NC(=O)C(C)(C)C)nc2NC(=O)C(C)(C)C)nc1": "25",
                "Brc1cc2cccnc2cc1": "26",
                "Brc1cnc(N2CCCCC2)nc1": "27",
                "Brc1ncc(c2ccccc2)s1": "28",
                "Brc1cscc1C": "29",
                "Brc1cccs1": "3",
                "Brc1c2c(cccc2)ncc1": "30",
                "Brc1cn(C(c2ccccc2)(c2ccccc2)c2ccccc2)nc1": "31",
                "Brc1c(C)cc(OC)cc1": "32",
                "Brc1ccc2c(c3ccccc3n2CC)c1": "33",
                "Brc1cn(cn1)C": "34",
                "Brc1cnc2n1cccc2": "35",
                "Brc1cnc2n1nccc2": "36",
                "Brc1cnc2n1ccnc2": "37",
                "Brc1csc(Oc2ccccc2)n1": "38",
                "Brc1ccccc1n1nccc1": "39",
                "Brc1ccc(C(F)(F)F)cc1": "4",
                "Brc1cn(c2ccccc2)nc1": "40",
                "Brc1nc(OC)ccc1": "41",
                "Brc1ccnc2c1ccn2S(=O)(=O)c1ccc(C)cc1": "42",
                "Brc1ccnc2c1cc[nH]2": "43",
                "Brc1cnc(c2ccccc2)o1": "44",
                "Brc1ncccc1": "45",
                "Brc1ccncc1": "46",
                "Brc1cc(C(F)(F)F)cc(C(F)(F)F)c1": "47",
                "Brc1c(C(F)(F)F)cccc1": "48",
                "Brc1ccc2c(cccc2)c1": "49",
                "Brc1cc2nccnc2cc1": "5",
                "Cn1cc(Br)ccc1=O": "50",
                "Brc1ccc(Oc2ccccc2)cc1": "6",
                "Brc1c(OC)cccc1": "7",
                "Brc1c2c(cccc2)ccc1": "8",
                "Brc1c(C)n(C)nc1C": "9",
                "c1cc(cc(c1)Br)C(=O)N(C)C": "2001",
                "c1cccc(c1CO[C@H]1CCCCO1)Br": "2002",
                "c12c(cc(cc1)Br)sc(n2)C": "2003",
                "C1COC(O1)c1cc(ccc1)Br": "2004",
                "c1(cccc2c1cnn2C)Br": "2005",
                "n1ccc(c2c1n(cc2)COCC[Si](C)(C)C)Br": "2006",
                "n1ccc(c2c1n(cc2)COC)Br": "2007",
                "n1ccc(c2c1n(cc2)Cc1ccccc1)Br": "2008",
                "O1C(=O)C(=C(C1(C)C)c1ccc(cc1)S(=O)(=O)C)Oc1ccc(cn1)Br": "2009",
                "n12-c3c(C(=O)N(Cc1c(nc2)C(=O)OCC)C)cc(cc3)Br": "2010",
                "Clc1cc2CCc3cc(cnc3C(=C3CCN(C(=O)OCC)CC3)c2cc1)Br": "2011",
                "c1cc(cc2c1n(c1c2cc(cc1)[C@H]1N(CCOC1)COC)CC)Br": "2012",
                "c1cc(cc2c1n(c1c2cc(cc1)[C@H]1N(CCOC1)S(=O)(=O)c1ccc(cc1)C)CC)Br": "2013",
                "c1cc(cc2c1n(c1c2cc(cc1)[C@H]1N(CCOC1)C(=O)OC(C)(C)C)CC)Br": "2014",
                "c1cc(cc2c1n(c1c2cc(cc1)[C@H]1N(CCOC1)Cc1ccc(cc1)OC)CC)Br": "2015"
            },
            'nuc': {
                "Nc1ccccc1": "1",
                "NC(=O)Cc1ccccc1": "10",
                "Nc1ccccn1": "11",
                "Cn1cc(N)cn1": "12",
                "CNC(=O)C": "13",
                "Cn1nccc1N": "14",
                "CNc1ccccn1": "15",
                "Cn1ccc(N)n1": "16",
                "c1cc[nH]c1": "17",
                "c1ccc2[nH]ccc2c1": "18",
                "Nc1cccnc1": "19",
                "Cc1cccc(C)c1N": "2",
                "c1c[nH]cn1": "20",
                "c1cn[nH]c1": "21",
                "CNc1cccnc1": "22",
                "N=C(c1ccccc1)c1ccccc1": "23",
                "NN=C(c1ccccc1)c1ccccc1": "24",
                "Nc1cnco1": "25",
                "Nc1ncco1": "26",
                "Nc1cncs1": "27",
                "Nc1cscn1": "28",
                "COc1cnc(N)s1": "29",
                "COc1ccc(N)cc1": "3",
                "Nc1ncc(C(F)(F)F)s1": "30",
                "NCC1CCCCC1": "31",
                "CN1CCC(CN)CC1": "32",
                "C[C@H](N)c1ccccc1": "34",
                "C[C@H]1CC[C@@H](N)CC1": "35",
                "C[C@H]1CC[C@H](N)CC1": "36",
                "CC(C)(C)N": "37",
                "C1CCNCC1": "38",
                "C[C@H]1CCC[C@H](C)N1": "39",
                "COc1cc(C)c(N)c(C)c1": "4",
                "C[C@@H]1CCC[C@H](C)N1": "40",
                "CC1(C)CCCC(C)(C)N1": "41",
                "C1C[C@@H]2CC[C@H](C1)N2": "42",
                "C1COCCN1": "43",
                "C1CCNC1": "44",
                "C1CC[C@@H]2CNCC[C@@H]2C1": "45",
                "C1CC[C@@H]2CNCC[C@H]2C1": "46",
                "CCCCNCCCC": "47",
                "CNCc1ccccc1": "48",
                "c1ccc(CNCc2ccccc2)cc1": "49",
                "Cc1ccc(N)cc1": "5",
                "CNC(C)(C)C": "50",
                "CNc1ccccc1": "6",
                "CCCC(=O)N": "7",
                "NC(=O)c1ccccc1": "8",
                "CC(C)(C)C(=O)N": "9",
                "O=C1NCCO1": "1001",
                "C1CCCNCC1": "1002",
                "Nc1cnccn1": "1003",
                "CC(C)(C)OC(=O)N": "1004",
                "C1CSCCN1": "1005",
                "CN1CCNCC1": "1006",
                "CC(C)(C)OC(=O)[C@@H]1CCC(=O)N1": "1007",
                "COc1ccc(CO[C@@H]2[C@@H](N)[C@H](Oc3ccc(OC)cc3)O[C@@H]3CO[C@@H](c4ccccc4)O[C@@H]23)cc1": "1008",
                "CC1(C)O[C@@H]2[C@@H](N)[C@@H]3OC[C@@H](O3)[C@@H]2O1": "1009",
                "Cc1c[nH]c2ccc(F)cc12": "1010",
                "COc1ccc(CO[C@@H](C)[C@H](N)c2cc3ccccc3s2)cc1": "1011"
            }
        }
        for config in configs:
            if config[mol_type] in ref[mol_type]:
                name_key = f'{mol_type}_name'
                name_map[config[name_key]] = ref[mol_type][config[mol_type]]
                config[name_key] = ref[mol_type][config[mol_type]]
                    
        return (configs, name_map)