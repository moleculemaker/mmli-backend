import io
import os
import json

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
    def validate_and_update_config(job_config: dict):
        reaction_sites, el_mol_str, _, _ = SomnService.check_user_input_substrates(job_config['el'], job_config['el_input_type'], 'el')
        
        if len(reaction_sites) > 1 or \
            job_config['el_input_type'] == 'cdxml' or \
            job_config['el_input_type'] == 'cml':
                job_config['el'] = el_mol_str
            
        reaction_sites, nuc_mol_str, _, _ = SomnService.check_user_input_substrates(job_config['nuc'], job_config['nuc_input_type'], 'nuc')
        
        if len(reaction_sites) > 1 or \
            job_config['nuc_input_type'] == 'cdxml' or \
            job_config['nuc_input_type'] == 'cml':
                job_config['nuc'] = nuc_mol_str
            
        return job_config
    
    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        job = await db.get(Job, job_id)
        job_info = json.loads(job.job_info)
        file_name = f"{job_info['nuc_name']}_{job_info['el_name']}_processed.csv"
        
        csv_content = service.get_file(bucket_name, f"/{job_id}/out/{file_name}")
        if csv_content is None:
            filename = "results/" + job_id + "/" + job_id + ".csv"
            raise HTTPException(status_code=404, detail=f"File {filename} not found")
        
        df = pd.read_csv(io.BytesIO(csv_content))
        retVal = []
        for index, row in df.iterrows():
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
            retVal.append(data)
            
        return retVal
    
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