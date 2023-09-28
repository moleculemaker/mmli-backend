from typing import List
from fastapi import HTTPException
import requests
import io
import os
import csv
import zipfile
import pandas as pd
from io import StringIO

from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.pubchem_service import PubChemService
from models.molecule import Molecule
from fastapi import Depends
from fastapi.responses import JSONResponse

class ChemScraperService:
    chemscraper_api_baseURL = os.environ.get("CHEMSCRAPER_API_BASE_URL")

    def __init__(self) -> None:
        pass

    def fetchExternalDataAndStoreResults(self, bucket_name: str, jobId: str, tsv_content: bytes,  service: MinIOService):
        reader = csv.reader(tsv_content.decode().splitlines(), delimiter='\t')

        doc_no = file_path = page_no = SMILE = minX = minY = maxX = maxY = SVG = PubChemCID = chemicalSafety = Description = None
        name = molecularFormula = molecularWeight = None
        molecules = []
        id = 0
        otherInstancesDict = {}
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
                    SVG = RDKitService.renderSVGFromSMILE(SMILE)
                    PubChemCID, name, molecularFormula, molecularWeight =  PubChemService.queryMoleculeProperties(SMILE)
                    location = " | page: " + page_no
                    if PubChemCID != 'Unavailable' and PubChemCID != '0':
                        chemicalSafety, Description = PubChemService.getAdditionalProperties(PubChemCID)
                    else:
                        chemicalSafety, Description = [], 'Unavailable'
                    if SMILE in otherInstancesDict:
                        otherInstancesDict[SMILE].append(page_no)
                    else:
                        otherInstancesDict[SMILE] = [page_no]
                    molecules.append(
                        Molecule(
                            id=id,
                            doc_no=doc_no,
                            file_path=file_path,
                            page_no=page_no,
                            name = name,
                            SMILE=SMILE, 
                            structure=SVG, 
                            minX=minX, 
                            minY=minY, 
                            width=maxX-minX, 
                            height=maxY-minY,
                            PubChemCID = PubChemCID,
                            molecularFormula = molecularFormula,
                            molecularWeight = molecularWeight,
                            chemicalSafety = chemicalSafety,
                            Description = Description,
                            Location = location,
                            OtherInstances = []
                        )
                    )
                    id += 1
        for molecule in molecules:
            pages = otherInstancesDict.get(molecule.SMILE, [])
            # pages = ', '.join(pages)
            # otherInstances = " | page(s): " + pages
            molecule.OtherInstances = pages

        data = [m.dict() for m in molecules]
        for d in data:
            d['chemicalSafety'] = ', '.join(d['chemicalSafety'])
            d['OtherInstances'] = ', '.join(d['OtherInstances'])
        df = pd.DataFrame(data)
        csv_buffer = StringIO()
        df.to_csv(csv_buffer)
        csv_data = csv_buffer.getvalue().encode('utf-8')

        upload_result = service.upload_file(bucket_name, "results/" + jobId + "/" + jobId + ".csv", csv_data)
        if not upload_result:
            raise HTTPException(status_code=500, detail="Unable to store CSV")
        return

    def runChemscraperOnDocument(self, bucket_name: str, filename: str, objectPath: str, jobId: str, service: MinIOService):
        data = service.get_file(bucket_name, objectPath)
        if data is None:
            raise HTTPException(status_code=404, detail="File not found")
        
        data_bytes = io.BytesIO(data)
        response = requests.post(self.chemscraper_api_baseURL + '/extractPdf', files={'pdf': (filename, data_bytes)})

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
                        self.fetchExternalDataAndStoreResults(bucket_name, jobId, tsv_data, service)
                        return True
        else:
            error_content = response.text.encode()
            upload_result = service.upload_file(bucket_name, "errors/" + jobId + ".txt", error_content)
        return False
