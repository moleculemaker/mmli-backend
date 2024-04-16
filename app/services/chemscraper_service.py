from typing import List
from fastapi import HTTPException
import requests
import io
import os
import csv
import zipfile
import time
import pandas as pd
from io import StringIO
from models.sqlmodel.models import Job
from models.enums import JobStatus

from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.pubchem_service import PubChemService
from services.email_service import EmailService
from models.molecule import Molecule
from fastapi import Depends
from fastapi.responses import JSONResponse

from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from models.sqlmodel.models import FlaggedMolecule

class ChemScraperService:
    chemscraper_api_baseURL = os.environ.get("CHEMSCRAPER_API_BASE_URL")
    chemscraper_frontend_baseURL = os.environ.get("CHEMSCRAPER_FRONTEND_URL")

    def __init__(self, db) -> None:
        self.db = db

    async def fetchExternalDataAndStoreResults(self, bucket_name: str, jobId: str, tsv_content: bytes,  service: MinIOService):
        reader = csv.reader(tsv_content.decode().splitlines(), delimiter='\t')
        rdkitService = RDKitService()
        doc_no = file_path = page_no = SMILE = minX = minY = maxX = maxY = SVG = None
        molecules = []
        id = 0
        otherInstancesDict = {}
        SMILE_LIST = []

        for row in reader:
            if not row:
                continue
            if row[0] == "D":
                doc_no, file_path = row[1], row[2]

            if row[0] == "P":
                page_no = row[1]

            if row[0] == "SMI":
                SMILE = row[2]
                minX, minY, maxX, maxY = map(int, row[3:7])
                if all([doc_no, file_path, page_no, SMILE, minX, minY, maxX, maxY]):
                    # Only molecules having all these fields available are processed
                    SMILE_LIST.append(SMILE)
                    svg_filename = f"Page_{page_no.zfill(3)}_No{row[1].zfill(3)}.svg"
                    SVG = service.get_file(bucket_name, "results/" + jobId + '/molecules/' + svg_filename)
                    if SVG is None:
                        print("SVG not found, generating using rdkit")
                        SVG = RDKitService.renderSVGFromSMILE(SMILE)

                    location = " | page: " + page_no

                    if SMILE in otherInstancesDict:
                        otherInstancesDict[SMILE].append(page_no)
                    else:
                        otherInstancesDict[SMILE] = [page_no]
                    try:
                        fingerprint = rdkitService.getFingerprint(SMILE)
                    except Exception as e:
                        print("Could not generate fingerprint for: " + SMILE)
                        fingerprint = "0"
                    try:
                        atom_count = rdkitService.getAtomCount(SMILE)
                    except Exception as e:
                        print("Could not generate fingerprint for: " + SMILE)
                        atom_count = 0
                    molecules.append(
                        Molecule(
                            id=id,
                            flagged=False,
                            atom_count=atom_count,
                            doc_no=doc_no,
                            file_path=file_path,
                            page_no=page_no,
                            SMILE=SMILE, 
                            structure=SVG, 
                            minX=minX, 
                            minY=minY, 
                            width=maxX-minX, 
                            height=maxY-minY,
                            Location = location,
                            OtherInstances = [],
                            fingerprint = fingerprint
                        )
                    )
                    id += 1

        # Only for debugging
        # TODO: Remove after Pub Chem Batching tested in PROD
        print('=== Printing Smile List ====')
        print(SMILE_LIST)
        print('=== End Printing Smile List ====')

        # Get data for all molecules
        pubChemService = PubChemService()
        molecules_data = await pubChemService.getDataForAllMolecules(SMILE_LIST)

        # Only for debugging
        # TODO: Remove after Pub Chem Batching tested in PROD
        print('======== Printing All Molecule Data =======')
        data_idx = 0
        while data_idx < len(molecules_data):
            print(molecules_data[data_idx], ' ', molecules_data[data_idx+1], ' ', molecules_data[data_idx+2], ' ',molecules_data[data_idx+3], ' ', molecules_data[data_idx+4])
            data_idx += 5
        print('======== End Printing All Molecule Data ======')

        # To iterate Molecule Data Array - molecules_data
        molecules_data_idx = 0

        # Setting Pubchem results directly to CSV
        data = [m.dict() for m in molecules]
        for d in data:
            d['chemicalSafety'] = ', '.join(d['chemicalSafety'])
            d['OtherInstances'] = ', '.join(otherInstancesDict.get(d['SMILE'], []))
            d['PubChemCID'] = molecules_data[molecules_data_idx + 1]
            d['molecularFormula'] = molecules_data[molecules_data_idx + 2]
            d['molecularWeight'] = molecules_data[molecules_data_idx + 3]
            d['name'] = molecules_data[molecules_data_idx + 4]
            molecules_data_idx += 5

        df = pd.DataFrame(data)
        csv_buffer = StringIO()
        df.to_csv(csv_buffer)
        csv_data = csv_buffer.getvalue().encode('utf-8')

        upload_result = service.upload_file(bucket_name, "results/" + jobId + "/" + jobId + ".csv", csv_data)
        if not upload_result:
            raise HTTPException(status_code=500, detail="Unable to store CSV")
        return

    async def update_job_phase(self, jobObject, phase: JobStatus):
        jobObject.phase = phase
        if phase == JobStatus.PROCESSING:
            jobObject.time_start = int(time.time())
        else:
            jobObject.time_end = int(time.time())
        self.db.add(jobObject)
        await self.db.commit()

    async def runChemscraperOnDocument(self, bucket_name: str, filename: str, objectPath: str, jobId: str, service: MinIOService, email_service: EmailService):
        # Get Job Object
        db_job : Job = await self.db.get(Job, jobId)

        # Update Job Status to Processing
        await self.update_job_phase(db_job, JobStatus.PROCESSING)

        data = service.get_file(bucket_name, objectPath)
        if data is None:
            await self.update_job_phase(db_job, JobStatus.ERROR)
            raise HTTPException(status_code=404, detail="File not found")
        
        data_bytes = io.BytesIO(data)
        params = {"generate_svg": True}
        response = requests.post(self.chemscraper_api_baseURL + '/extractPdf', files={'pdf': (filename, data_bytes)}, params=params)

        if response.status_code == 200:
            zip_file = zipfile.ZipFile(io.BytesIO(response.content))
            file_list = zip_file.namelist()
            tsv_file_name = next((file for file in file_list if file.endswith('.tsv')), None)

            for filename in file_list:
                file_info = zip_file.getinfo(filename)
                if file_info.is_dir():
                    continue
                if not filename.endswith('.tsv'):
                    file_data = zip_file.read(filename)
                    upload_result = service.upload_file(bucket_name, "results/" + jobId + '/' + filename, file_data)

            if tsv_file_name is not None:
                with zip_file.open(tsv_file_name) as tsv_file:
                    tsv_data = tsv_file.read()
                    upload_result = service.upload_file(bucket_name, "results/" + jobId + '/' + jobId + ".tsv", tsv_data)
                    if upload_result:
                        await self.fetchExternalDataAndStoreResults(bucket_name, jobId, tsv_data, service)
                        await self.update_job_phase(db_job, JobStatus.COMPLETED)
                        if(db_job.email):
                            try:
                                email_service.send_email(db_job.email, f'''Result for your ChemScraper Job ({db_job.job_id}) is ready''', f'''The result for your ChemScraper Job is available at {self.chemscraper_frontend_baseURL}/results/{db_job.job_id}''')
                            except Exception as e:
                                print(e)
                        return True
        else:
            error_content = response.text.encode()
            upload_result = service.upload_file(bucket_name, "errors/" + jobId + ".txt", error_content)
            await self.update_job_phase(db_job, JobStatus.ERROR)
            if(db_job.email):
                try:
                    email_service.send_email(db_job.email, f'''ChemScraper Job ({db_job.job_id}) failed''', f'''An error occurred in computing the result for your ChemScraper job.''')
                except Exception as e:
                    print(e)
        return False

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        rdkitService = RDKitService()
        csv_content = service.get_file(bucket_name, "results/" + job_id + "/" + job_id + ".csv")
        if csv_content is None:
            filename = "results/" + job_id + "/" + job_id + ".csv"
            raise HTTPException(status_code=404, detail=f"File {filename} not found")
        molecules = []
        df = pd.read_csv(io.BytesIO(csv_content))

        # Returns True is smile is flagged
        def is_row_flagged(job_id, row, flagged_molecules):
            for mol in flagged_molecules.scalars().all():
                if row['SMILE'] == mol.smile:
                    return True
            return False

        for index, row in df.iterrows():
            smile = row['SMILE']
            doc_id = row['doc_no']

            flagged_molecules = await db.execute(select(FlaggedMolecule).where(
                FlaggedMolecule.smile == smile
                and FlaggedMolecule.job_id == job_id
                and FlaggedMolecule.doc_id == doc_id))

            # Convert the 'chemicalSafety' and 'OtherInstances' strings back into lists
            chemicalSafety = str(row['chemicalSafety']).split(', ')
            OtherInstances = str(row['OtherInstances']).split(', ')

            # Create a Molecule object and append it to the list
            molecule = Molecule(id=row['id'],
                                doc_no=row['doc_no'],
                                atom_count=row['atom_count'] if 'atom_count' in row else rdkitService.getAtomCount(smile),
                                file_path=row['file_path'],
                                page_no=row['page_no'],
                                name=row['name'],
                                flagged=is_row_flagged(job_id, row, flagged_molecules),
                                SMILE=row['SMILE'],
                                structure=row['structure'],
                                minX=row['minX'],
                                minY=row['minY'],
                                width=row['width'],
                                height=row['height'],
                                PubChemCID=row['PubChemCID'],
                                molecularFormula=row['molecularFormula'],
                                molecularWeight=row['molecularWeight'],
                                chemicalSafety=chemicalSafety,
                                Description=row['Description'],
                                Location=row['Location'],
                                OtherInstances=OtherInstances,
                                fingerprint=row['fingerprint'])
            molecules.append(molecule)
        return molecules