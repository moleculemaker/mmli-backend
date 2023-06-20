import csv
import uuid
from typing import List
from datetime import datetime
from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, status, Query, BackgroundTasks
from fastapi.responses import JSONResponse

from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.pubchem_service import PubChemService
from services.chemscraper_service import ChemScraperService

from models.molecule import Molecule
from models.analyzeRequestBody import AnalyzeRequestBody
from typing import Optional

router = APIRouter()

@router.post("/{bucket_name}/upload")
async def upload_file(bucket_name: str, file: UploadFile = File(...), job_id: Optional[str] = Query(None), service: MinIOService = Depends()):
    first_four_bytes = file.file.read(4)
    file.file.seek(0)
    if first_four_bytes == b'%PDF':
        if job_id == "":
            job_id = str(uuid.uuid4())
        file_content = await file.read()
        upload_result = service.upload_file(bucket_name, "inputs/" + job_id + '/' + file.filename, file_content)
        if upload_result:
            content = {"jobID": job_id, "uploaded_at": datetime.now().isoformat()}
            return JSONResponse(content=content, status_code=status.HTTP_200_OK)    
    return JSONResponse(content={"error": "Unable to upload file"}, status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)

@router.post("/chemscraper/analyze")
async def analyze_documents(requestBody: AnalyzeRequestBody, background_tasks: BackgroundTasks, service: MinIOService = Depends()):
    # Analyze only one document for NSF demo
    if len(requestBody.fileList) > 0 and requestBody.jobId != "":
        filename = requestBody.fileList[0]
        chemscraperService = ChemScraperService()
        objectPath = f"inputs/{requestBody.jobId}/{filename}"
        background_tasks.add_task(chemscraperService.runChemscraperOnDocument, 'chemscraper', filename, objectPath, requestBody.jobId, service)
        # chemscraperService.runChemscraperOnDocument('chemscraper', filename, objectPath, requestBody.jobId, service)
        content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
        return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 

@router.get("/{bucket_name}/result-status/{job_id}")
def get_result_status(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    result_status = service.check_file_exists(bucket_name, "results/" + job_id + ".tsv")
    error_status = service.check_file_exists(bucket_name, "errors/" + job_id + ".txt")
    if result_status:
        # JSONResponse(content={"Result": "Ready"}, status_code=status.HTTP_200_OK)
        return "Ready"
    elif error_status:
        # JSONResponse(content={"Result": "Error"}, status_code=status.HTTP_200_OK)
        return "Error"
    else:
        # JSONResponse(content={"Result": "Processing"}, status_code=status.HTTP_200_OK)
        return "Result"

@router.get("/{bucket_name}/results/{job_id}",response_model=List[Molecule])
def get_results(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    tsv_content = service.get_file(bucket_name, "results/" + job_id + ".tsv")
    if tsv_content is None:
        raise HTTPException(status_code=404, detail="File not found") 
    
    reader = csv.reader(tsv_content.decode().splitlines(), delimiter='\t')

    doc_no = file_path = page_no = SMILE = minX = minY = maxX = maxY = SVG = PubChemCID = chemicalSafety = Description = None
    name = molecularFormula = molecularWeight = None
    molecules = []
    id = 0
    for row in reader:
        if not row:
            continue
        if row[0] == "D":
            doc_no, file_path = int(row[1]), row[2]

        if row[0] == "P":
            page_no = int(row[1])

        if row[0] == "SMI":
            SMILE = row[2]
            minX, minY, maxX, maxY = map(int, row[3:7])
            if all([doc_no, file_path, page_no, SMILE, minX, minY, maxX, maxY]):
                SVG = RDKitService.renderSVGFromSMILE(SMILE)
                PubChemCID, name, molecularFormula, molecularWeight =  PubChemService.queryMoleculeProperties(SMILE)
                if PubChemCID != 'Unknown' and PubChemCID != '0':
                    chemicalSafety, Description = PubChemService.getAdditionalProperties(PubChemCID)
                else:
                    chemicalSafety, Description = [], 'Unknown'
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
                        Description = Description
                    )
                )
                id += 1
    return molecules

@router.get("/{bucket_name}/inputs/{job_id}")
def get_input_file(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    pdf_urls = service.get_file_urls(bucket_name, "inputs/" + job_id + "/")
    if pdf_urls is None:
        raise HTTPException(status_code=404, detail="File not found") 
    return pdf_urls

@router.get("/{bucket_name}/errors/{job_id}")
def get_errors(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    error_content = service.get_file(bucket_name, "errors/" + job_id + ".txt")
    if error_content is None:
        raise HTTPException(status_code=404, detail="File not found") 
    return error_content
