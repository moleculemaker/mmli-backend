import csv
import io

import uuid
import zipfile
from typing import List
from datetime import datetime

import pandas as pd
from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, status
from fastapi.responses import JSONResponse
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession
from starlette.responses import FileResponse

from models.exportRequestBody import ExportRequestBody
from models.sqlmodel.db import get_session
from models.sqlmodel.models import FlaggedMolecule
from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.pubchem_service import PubChemService

from models.molecule import Molecule
from typing import Optional

router = APIRouter()


@router.post("/{bucket_name}/upload", tags=['Files'])
async def upload_file(bucket_name: str, file: UploadFile = File(...), job_id: Optional[str] = "", minio: MinIOService = Depends()):
    first_four_bytes = file.file.read(4)
    file.file.seek(0)
    if first_four_bytes == b'%PDF':
        if job_id == "":
            job_id = str(uuid.uuid4())

        file_content = await file.read()
        upload_result = minio.upload_file(bucket_name, "inputs/" + job_id + '/' + file.filename, file_content)
        if upload_result:
            content = {"jobID": job_id, "uploaded_at": datetime.now().isoformat()}
            return JSONResponse(content=content, status_code=status.HTTP_200_OK)
    return JSONResponse(content={"error": "Unable to upload file"}, status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)


# Returns True is smile is flagged
def is_row_flagged(job_id, row, flagged_molecules):
    for mol in flagged_molecules.scalars().all():
        if row['SMILE'] == mol.smile:
            return True
    return False


@router.get("/{bucket_name}/results/{job_id}", response_model=List[Molecule], tags=['Files'])
async def get_results(bucket_name: str, job_id: str, service: MinIOService = Depends(), db: AsyncSession = Depends(get_session)):
    csv_content = service.get_file(bucket_name, "results/" + job_id + "/" + job_id + ".csv")
    if csv_content is None:
        filename = "results/" + job_id + "/" + job_id + ".csv"
        raise HTTPException(status_code=404, detail=f"File {filename} not found")
    molecules = []
    df = pd.read_csv(io.BytesIO(csv_content))

    for index, row in df.iterrows():
        doc_id = row['doc_no']

        flagged_molecules = await db.execute(select(FlaggedMolecule).where(
            FlaggedMolecule.job_id == job_id
            and FlaggedMolecule.doc_id == doc_id))

        # Convert the 'chemicalSafety' and 'OtherInstances' strings back into lists
        chemicalSafety = str(row['chemicalSafety']).split(', ')
        OtherInstances = str(row['OtherInstances']).split(', ')

        # Create a Molecule object and append it to the list
        molecule = Molecule(id=row['id'],
                            doc_no=row['doc_no'],
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


@router.get("/{bucket_name}/inputs/{job_id}", tags=['Files'])
def get_input_file(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    pdf_urls = service.get_file_urls(bucket_name, "inputs/" + job_id + "/")
    if pdf_urls is None:
        filename = "inputs/" + job_id + "/"
        raise HTTPException(status_code=404, detail=f"Files {filename} not found")
    return pdf_urls


@router.get("/{bucket_name}/errors/{job_id}", tags=['Files'])
def get_errors(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    error_content = service.get_file(bucket_name, "errors/" + job_id + ".txt")
    if error_content is None:
        filename = "errors/" + job_id + ".txt"
        raise HTTPException(status_code=404, detail=f"File {filename} not found")
    return error_content


@router.post("/{bucket_name}/export-results", tags=['Files'])
async def analyze_documents(bucket_name: str, requestBody: ExportRequestBody, service: MinIOService = Depends()):
    # Analyze only one document for NSF demo
    if requestBody.jobId == "":
        raise HTTPException(status_code=404, detail="Invalid Job ID")
    if requestBody.jobId != "":
        objectPathPrefix = "results/" + requestBody.jobId + "/"
        files_count = 0
        filename = f'chemscraper_{requestBody.jobId}.zip'
        with zipfile.ZipFile(filename, "w") as new_zip:
            if requestBody.cdxml:
                if requestBody.cdxml_filter == "all_molecules":
                    object_path = objectPathPrefix + "molecules_full_cdxml/molecules.cdxml"
                    cdxml_file_data = service.get_file(bucket_name, object_path)
                    if cdxml_file_data is None:
                        filename = objectPathPrefix + "molecules_full_cdxml/molecules.cdxml"
                        raise HTTPException(status_code=404, detail=f"File {filename} not found")
                    new_zip.writestr(requestBody.jobId + ".cdxml", cdxml_file_data)
                    files_count += 1
                elif requestBody.cdxml_filter == "single_page" and len(requestBody.cdxml_selected_pages) > 0:
                    cdxml_file_data = service.get_file(bucket_name,
                                                       objectPathPrefix + "molecules_all_pages/Page_" + f"{requestBody.cdxml_selected_pages[0]:03d}" + ".cdxml")
                    if cdxml_file_data is None:
                        filename = objectPathPrefix + "molecules_all_pages/Page_" + f"{requestBody.cdxml_selected_pages[0]:03d}" + ".cdxml"
                        raise HTTPException(status_code=404, detail=f"File {filename} not found")
                    new_zip.writestr(requestBody.jobId + "_Page_" + f"{requestBody.cdxml_selected_pages[0]:03d}" + ".cdxml",
                                     cdxml_file_data)
                    files_count += 1
            if requestBody.csv:
                if requestBody.csv_filter == "full_table":
                    csv_file_data = service.get_file(bucket_name, objectPathPrefix + requestBody.jobId + ".csv")
                    if csv_file_data is None:
                        filename = objectPathPrefix + requestBody.jobId + ".csv"
                        raise HTTPException(status_code=404, detail=f"File {filename} not found")
                    new_zip.writestr(requestBody.jobId + ".csv", csv_file_data)
                    files_count += 1
                elif requestBody.csv_filter == "current_view":
                    csv_file_data = service.get_file(bucket_name, objectPathPrefix + requestBody.jobId + ".csv")
                    if csv_file_data is None:
                        filename = objectPathPrefix + requestBody.jobId + ".csv"
                        raise HTTPException(status_code=404, detail=f"File {filename} not found")
                    csvfile = io.StringIO(csv_file_data.decode('utf-8'))
                    reader = csv.DictReader(csvfile)
                    rows = [row for row in reader]
                    reordered_rows = [rows[i] for i in requestBody.csv_molecules]
                    output_csv = io.StringIO()
                    writer = csv.DictWriter(output_csv, fieldnames=reordered_rows[0].keys())
                    writer.writeheader()
                    writer.writerows(reordered_rows)
                    new_zip.writestr(requestBody.jobId + ".csv", output_csv.getvalue())
                    files_count += 1

        if files_count > 0:
            return FileResponse(filename, media_type='application/zip', filename=filename)
        else:
            raise HTTPException(status_code=400, detail="Bad Request")
