import csv
import io

import uuid
import zipfile
from datetime import datetime

from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, status
from fastapi.responses import JSONResponse
from sqlmodel.ext.asyncio.session import AsyncSession
from starlette.responses import FileResponse

from config import get_logger
from models.exportRequestBody import ExportRequestBody
from models.sqlmodel.db import get_session

from models.enums import JobType
from services.molli_service import MolliService
from services.clean_service import CleanService

from services.novostoic_service import NovostoicService
from services.somn_service import SomnService
from services.minio_service import MinIOService
from services.chemscraper_service import ChemScraperService


from typing import Optional

router = APIRouter()

log = get_logger(__name__)

@router.post("/{bucket_name}/upload", tags=['Files'])
async def upload_file(bucket_name: str, file: UploadFile = File(...), job_id: Optional[str] = "", minio: MinIOService = Depends()):
    first_four_bytes = file.file.read(4)
    file.file.seek(0)
    if bucket_name != 'chemscraper' or first_four_bytes == b'%PDF':
        if job_id == "":
            job_id = str(uuid.uuid4()).replace('-', '')

        file_content = await file.read()
        upload_result = minio.upload_file(bucket_name, job_id + '/in/' + file.filename, file_content)
        if upload_result:
            content = {"jobID": job_id, "uploaded_at": datetime.now().isoformat()}
            return JSONResponse(content=content, status_code=status.HTTP_200_OK)

    log.error(f'Failed to upload file: {file.filename}')
    return JSONResponse(content={"error": "Unable to upload file"}, status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)


@router.get("/{bucket_name}/results/{job_id}", tags=['Files'])
async def get_results(bucket_name: str, job_id: str, service: MinIOService = Depends(), db: AsyncSession = Depends(get_session)):
    if bucket_name == JobType.CHEMSCRAPER:
        print("Getting CHEMSCRAPER job result")
        return await ChemScraperService.resultPostProcess(bucket_name, job_id, service, db)
    
    elif bucket_name == JobType.MOLLI:
        print("Getting MOLLI job result")
        return await MolliService.molliResultPostProcess(bucket_name, job_id, service, db)

    elif bucket_name == JobType.CLEAN:
        print("Getting CLEAN job result")
        return await CleanService.cleanResultPostProcess(bucket_name, job_id, service, db)
        
    elif bucket_name == JobType.NOVOSTOIC_OPTSTOIC:
        print("Getting novostoic-optstoic job result")
        return await NovostoicService.optstoicResultPostProcess(bucket_name, job_id, service, db)
    
    elif bucket_name == JobType.NOVOSTOIC_PATHWAYS:
        print("Getting novostoic-pathways job result")
        return await NovostoicService.novostoicResultPostProcess(bucket_name, job_id, service, db)
        
    elif bucket_name == JobType.NOVOSTOIC_ENZRANK:
        print("Getting novostoic-enzrank job result")
        return await NovostoicService.enzRankResultPostProcess(bucket_name, job_id, service, db)

    elif bucket_name == JobType.NOVOSTOIC_DGPREDICTOR:
        print("Getting novostoic-dgpredictor job result")
        return await NovostoicService.dgPredictorResultPostProcess(bucket_name, job_id, service, db)

    elif bucket_name == JobType.SOMN:
        return await SomnService.resultPostProcess(bucket_name, job_id, service, db)

    else:
        raise HTTPException(status_code=400, detail="Invalid job type: " + bucket_name)


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
