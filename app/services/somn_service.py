import io
import os
import json

from fastapi import HTTPException
import pandas as pd
from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.models import Job

from services.minio_service import MinIOService

from typing import List, Literal
from rdkit import Chem
from openbabel import pybel as pb

class SomnException(Exception):
    pass

class SomnService:
    somn_frontend_baseURL = os.environ.get("SOMN_FRONTEND_URL")
    
    @staticmethod
    def gen3d_test(
        user_input: str, 
        input_type: Literal['smi', 'cml', 'cdxml']
    ):
        try:            
            obmol = pb.readstring(input_type, user_input)
            obmol.addh()
            
            # this step may fail, so we know SOMN cannot compute on the input
            obmol.make3D() 
            
            return obmol.write("mol2")
            
        except Exception as e:
            raise SomnException(f"Unable to generate 3D coordinates for {user_input}")
    
    @staticmethod
    def validate_and_update_config(job_config: dict):
        el_mol_str = SomnService.gen3d_test(job_config['el'], job_config['el_input_type'])
        reaction_sites = SomnService.check_user_input_substrates(job_config['el'], job_config['el_input_type'], 'el')
        
        if len(reaction_sites) > 1 or \
            job_config['el_input_type'] == 'cdxml' or \
            job_config['el_input_type'] == 'cml':
                job_config['el'] = el_mol_str
            
        nuc_mol_str = SomnService.gen3d_test(job_config['nuc'], job_config['nuc_input_type'])
        reaction_sites = SomnService.check_user_input_substrates(job_config['nuc'], job_config['nuc_input_type'], 'nuc')
        
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
        Verifies user input substrate, and if verification is successful, returns reaction sites if multiple.
        If everything is "normal", i.e., the user doesn't need to tell us more information, then it returns 0.
        If NONE are found for halides or reactive nitrogens, then this returns None.
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

        # add gen3d_test here to prevent the user from 
        # submitting a molecule that cannot be processed by SOMN
        obmol = SomnService.gen3d_test(user_input, input_type)

        # generate rdkit mol using obmol because rdkit
        # might have issues generating mol from smiles directly
        mol = Chem.MolFromMol2Block(obmol, sanitize=False, removeHs=False)

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
            return indices

        elif role.startswith('el'):
            if len(bromides) + len(chlorides) == 0:
                raise SomnException("No Br or Cl sites detected in electrophile!")
            
            try:
                chl_idxes = check_halides_aromatic(mol,chlorides)
                brm_idxes = check_halides_aromatic(mol,bromides)

                return chl_idxes + brm_idxes
            
            except Exception as e:
                raise e

        return []