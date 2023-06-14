import csv
import uuid
from typing import List
from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, status
from fastapi.responses import JSONResponse
from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.pubchem_service import PubChemService
from models.molecule import Molecule

router = APIRouter()
bucket_name = "chemscraper"

@router.post("/submit/{bucket_name}")
async def upload_file(file: UploadFile = File(...), service: MinIOService = Depends()):

    first_four_bytes = file.file.read(4)
    file.file.seek(0)

    if first_four_bytes == b'%PDF':
        job_id = uuid.uuid4()
        file_content = await file.read()
        upload_result = service.upload_file(bucket_name, "input/" + job_id + ".pdf", file_content)
        if upload_result:
            content = {"job_id": job_id}
            return JSONResponse(content=content, status_code=status.HTTP_200_OK)    
    return JSONResponse(content={"error": "Unable to upload file"}, status_code=status.HTTP_500_INTERNAL_SERVER_ERROR)

@router.get("/results/{bucket_name}/{job_id}",response_model=List[Molecule])
def get_results(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    tsv_content = service.get_file(bucket_name, "results/" + job_id + ".tsv")
    if tsv_content is None:
        raise HTTPException(status_code=404, detail="File not found") 
    
    reader = csv.reader(tsv_content.decode().splitlines(), delimiter='\t')

    doc_no = file_path = page_no = SMILE = minX = minY = maxX = maxY = SVG = PubChemCID = chemicalSafety = Description = None
    name = molecularFormula = molecularWeight = None
    molecules = []

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
                break
    return molecules
