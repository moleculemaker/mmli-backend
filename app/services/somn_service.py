import io
import os
import time
import asyncio
import json

from fastapi import HTTPException
import pandas as pd
from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.models import Job
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.email_service import EmailService

from typing import List, Literal
from rdkit import Chem
from openbabel import pybel as pb
from openbabel import openbabel as ob

class SomnException(Exception):
    pass

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
    def ob_test(smiles: str):
        """
        Generate 3D coordinates for a molecule from its SMILES string using OpenBabel.
        
        Parameters:
            smiles (str): SMILES representation of the molecule
            
        Returns:
            tuple[str, bool]: A tuple containing:
                - The molecule in MOL2 format with 3D coordinates
                - Boolean indicating if the molecule has chiral centers
                
        Raises:
            SomnException: If 3D coordinate generation fails
        """
        try:
            obmol = pb.readstring("smi", smiles)
        except Exception as e:
            raise SomnException(f"Invalid SMILES string: {smiles}")
        
        obmol.addh()
        
        try:            
            # this step may fail, so we know SOMN cannot compute on the input
            obmol.make3D() 
            
            return obmol.write("mol2"), SomnService.has_chiral(obmol)
            
        except Exception as e:
            raise SomnException(f"Unable to generate 3D coordinates for {smiles}")
    
    @staticmethod
    def validate_and_update_config(job_config: dict):
        reaction_sites, el_mol_str, _ = SomnService.check_user_input_substrates(job_config['el'], 'el')
        
        if len(reaction_sites) > 1:
            job_config['el'] = el_mol_str
            
        reaction_sites, nuc_mol_str, _ = SomnService.check_user_input_substrates(job_config['nuc'], 'nuc')
        
        if len(reaction_sites) > 1:
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
    def check_user_input_substrates(user_input, role: str):
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
            
            raise Exception("No reactive nitrogens detected in nucleophile!")
            
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
        (mol2, has_chiral) = SomnService.ob_test(user_input)

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
                    raise SomnException("Bromides detected in nucleophile!")

            if len(chlorides) != 0:
                aromatic_halides = check_halides_aromatic(mol,bromides)
                if any(idx in aromatic_halides for idx in chlorides):
                    raise SomnException("Chlorides detected in nucleophile!")
            
            if len(nitrogens) == 0:
                raise SomnException("No nitrogens detected in nucleophile!")
            
            indices = get_amine_ref_ns(mol,nitrogens)
            return (indices, mol2, has_chiral)

        elif role.startswith('el'):
            if len(bromides) + len(chlorides) == 0:
                raise SomnException("No Br or Cl sites detected in electrophile!")
            
            try:
                chl_idxes = check_halides_aromatic(mol,chlorides)
                brm_idxes = check_halides_aromatic(mol,bromides)

                return (chl_idxes + brm_idxes, mol2, has_chiral)
            
            except Exception as e:
                raise e

        return ([], "", has_chiral)
