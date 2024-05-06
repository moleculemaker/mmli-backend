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

from typing import List
from rdkit import Chem

class SomnService:
    somn_frontend_baseURL = os.environ.get("SOMN_FRONTEND_URL")

    def __init__(self, db) -> None:
        self.db = db

    async def update_job_phase(self, jobObject, phase: JobStatus):
        jobObject.phase = phase
        if phase == JobStatus.PROCESSING:
            jobObject.time_start = int(time.time())
        else:
            jobObject.time_end = int(time.time())
        self.db.add(jobObject)
        await self.db.commit()

    async def runSomn(self, bucket_name: str, payload, service: MinIOService, email_service: EmailService):
        job_id = payload["job_id"]
        db_job = await self.db.get(Job, job_id)

        await self.update_job_phase(db_job, JobStatus.PROCESSING)

        #TODO: use somn image
        await asyncio.sleep(3)

        sample_response = "somn:pong"
        upload_result = service.upload_file(bucket_name, f"results/{job_id}/{job_id}.csv", bytes(sample_response, 'utf-8'))

        if upload_result:
            await self.update_job_phase(db_job, JobStatus.COMPLETED)
            if(db_job.email):
                try:
                    email_service.send_email(
                        db_job.email, 
                        f'''Result for your Somn Job ({db_job.job_id}) is ready''', 
                        f'''The result for your Somn Job is available at {self.somn_frontend_baseURL}/results/{db_job.job_id}''')
                except Exception as e:
                    print(e)
                return True
        else:
            error_content = 'error'
            upload_result = service.upload_file(bucket_name, f"results/{job_id}/error.txt", bytes(error_content, 'utf-8'))
            await self.update_job_phase(db_job, JobStatus.ERROR)
            if db_job.email:
                try:
                    email_service.send_email(
                        db_job.email, 
                        f'''Somn Job ({db_job.job_id}) failed''', 
                        f'''An error occurred in computing the result for your Somn job.''')
                except Exception as e:
                    print(e)
            
        return False
    
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
    def check_user_input_substrates(user_input, role):
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
            for idx in ref_atom_idxes:
                neighbors = mol.GetAtomWithIdx(idx).GetNeighbors()
                n_hs = [1 for nbr in neighbors if nbr.GetSymbol() == "H"]
                reactive_nitrogens.append(sum(n_hs))

            reactive_nitrogens = np.array(reactive_nitrogens)
            atm_indices = np.where((reactive_nitrogens == np.max(reactive_nitrogens)) & (reactive_nitrogens != 0))[0]

            if len(atm_indices) >= 1:
                return atm_indices
            
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

        mol = Chem.MolFromSmiles(user_input)
        mol = Chem.AddHs(mol)

        bromides = get_atoms_by_symbol(mol, symbol="Br")
        chlorides = get_atoms_by_symbol(mol, symbol="Cl")
        nitrogens = get_atoms_by_symbol(mol, symbol="N")

        if role == 'nucleophile':
            import numpy as np
            if len(bromides) != 0:
                aromatic_halides = check_halides_aromatic(mol,bromides)
                if any(idx in aromatic_halides for idx in bromides):
                    raise Exception("Bromides detected in nucleophile!")

            if len(chlorides) != 0:
                aromatic_halides = check_halides_aromatic(mol,bromides)
                if any(idx in aromatic_halides for idx in chlorides):
                    raise Exception("Chlorides detected in nucleophile!")
            
            if len(nitrogens) == 0:
                raise Exception("No nitrogens detected in nucleophile!")
            
            indices = get_amine_ref_ns(mol,nitrogens)
            return indices

        elif role == 'electrophile':
            if len(bromides) + len(chlorides) == 0:
                raise Exception("No Br or Cl sites detected in electrophile!")
            
            try:
                chl_idxes = check_halides_aromatic(mol,chlorides)
                brm_idxes = check_halides_aromatic(mol,bromides)

                if len(brm_idxes + chl_idxes) == 0: 
                    return None
                
                return chl_idxes + brm_idxes
            
            except Exception as e:
                raise e

        return None
