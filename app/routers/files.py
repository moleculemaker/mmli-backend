import csv

import uuid
from typing import List
from datetime import datetime
from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, status
from fastapi.responses import JSONResponse

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


@router.get("/{bucket_name}/result-status/{job_id}")
def get_result_status(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    result_status = service.check_file_exists(bucket_name, "results/" + job_id + '/' + job_id + ".csv")
    error_status = service.check_file_exists(bucket_name, "errors/" + job_id + ".txt")
    if result_status:
        # JSONResponse(content={"Result": "Ready"}, status_code=status.HTTP_200_OK)
        return "Ready"
    elif error_status:
        # JSONResponse(content={"Result": "Error"}, status_code=status.HTTP_200_OK)
        return "Error"
    else:
        # JSONResponse(content={"Result": "Processing"}, status_code=status.HTTP_200_OK)
        return "Processing"


@router.get("/{bucket_name}/results/{job_id}", response_model=List[Molecule])
def get_results(bucket_name: str, job_id: str, service: MinIOService = Depends()):
    csv_content = service.get_file(bucket_name, "results/" + job_id + "/" + job_id + ".csv")
    if csv_content is None:
        raise HTTPException(status_code=404, detail="File not found")
    molecules = []
    df = pd.read_csv(io.BytesIO(csv_content))

    for index, row in df.iterrows():
        # Convert the 'chemicalSafety' and 'OtherInstances' strings back into lists
        chemicalSafety = str(row['chemicalSafety']).split(', ')
        OtherInstances = str(row['OtherInstances']).split(', ')

        # Create a Molecule object and append it to the list
        molecule = Molecule(id=row['id'],
                            doc_no=row['doc_no'],
                            file_path=row['file_path'],
                            page_no=row['page_no'],
                            name=row['name'],
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
                            OtherInstances=OtherInstances)
        molecules.append(molecule)
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


@router.post("/{bucket_name}/export-results")
async def analyze_documents(bucket_name: str, requestBody: ExportRequestBody, service: MinIOService = Depends()):
    # Analyze only one document for NSF demo
    if requestBody.jobId == "":
        raise HTTPException(status_code=404, detail="Invalid Job ID")
    if requestBody.jobId != "":
        objectPathPrefix = "results/" + requestBody.jobId + "/"
        files_count = 0
        filename = f'chemscraper_{requestBody.jobId}.zip'
        with zipfile.ZipFile(filename, "w") as new_zip:
            if (requestBody.cdxml):
                if (requestBody.cdxml_filter == "all_molecules"):
                    cdxml_file_data = service.get_file(bucket_name,
                                                       objectPathPrefix + "molecules_full_cdxml/molecules_allpages.cdxml")
                    new_zip.writestr(requestBody.jobId + ".cdxml", cdxml_file_data)
                    files_count += 1
                elif (requestBody.cdxml_filter == "single_page" and len(requestBody.cdxml_selected_pages) > 0):
                    cdxml_file_data = service.get_file(bucket_name,
                                                       objectPathPrefix + "molecules_all_pages/Page_" + str(
                                                           requestBody.cdxml_selected_pages[0]) + "_all.cdxml")
                    new_zip.writestr(requestBody.jobId + ".cdxml", cdxml_file_data)
                    files_count += 1
            if (requestBody.csv):
                if (requestBody.csv_filter == "full_table"):
                    csv_file_data = service.get_file(bucket_name, objectPathPrefix + requestBody.jobId + ".csv")
                    new_zip.writestr(requestBody.jobId + ".csv", csv_file_data)
                    files_count += 1
                elif (requestBody.csv_filter == "current_view"):
                    csv_file_data = service.get_file(bucket_name, objectPathPrefix + requestBody.jobId + ".csv")
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

        if (files_count > 0):
            return FileResponse(filename, media_type='application/zip', filename=filename)
        else:
            raise HTTPException(status_code=400, detail="Bad Request")
