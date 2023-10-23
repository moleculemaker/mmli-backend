
from datetime import datetime
from fastapi import APIRouter, Depends, status, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse

from services.minio_service import MinIOService
from services.rdkit_service import RDKitService
from services.chemscraper_service import ChemScraperService
from services.hCaptcha_service import HCaptchaService

from models.analyzeRequestBody import AnalyzeRequestBody

router = APIRouter()


@router.post("/chemscraper/analyze", tags=['ChemScraper'])
async def analyze_documents(requestBody: AnalyzeRequestBody, background_tasks: BackgroundTasks, service: MinIOService = Depends()):
    # Analyze only one document for NSF demo
    if len(requestBody.fileList) > 0 and requestBody.jobId != "":
        hcaptchaService = HCaptchaService()
        if(hcaptchaService.verify_captcha(requestBody.captcha_token)):
            filename = requestBody.fileList[0]
            chemscraperService = ChemScraperService()
            objectPath = f"inputs/{requestBody.jobId}/{filename}"
            background_tasks.add_task(chemscraperService.runChemscraperOnDocument, 'chemscraper', filename, objectPath, requestBody.jobId, service)
            content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
            return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED)
        else:
            raise HTTPException(status_code=400, detail="Could not verify CAPTCHA")

@router.get("/chemscraper/similarity-sorted-order/{job_id}")
def get_similarity_sorted_order(job_id: str, smile_string: str, service: MinIOService = Depends()):
    bucket_name = "chemscraper"
    csv_content = service.get_file(bucket_name, "results/" + job_id + "/" + job_id + ".csv")
    if csv_content is None:
        filename = "results/" + job_id + "/" + job_id + ".csv"
        raise HTTPException(status_code=404, detail=f"File {filename} not found")
    df = pd.read_csv(io.BytesIO(csv_content))
    # Check if sort_column exists in DataFrame
    sort_column = "fingerprint"
    if sort_column not in df.columns:
        raise HTTPException(status_code=400, detail=f"Column {sort_column} not found in CSV file")
    rdkitService = RDKitService()
    input_fingerprint =  rdkitService.getFingerprint(smile_string)
    # Calculate Tanimoto similarity for each molecule in the DataFrame
    similarity_scores = []
    for index, row in df.iterrows():
        fingerprint = row['fingerprint']
        similarity = rdkitService.getTanimotoSimilarity(input_fingerprint, fingerprint)
        similarity_scores.append(similarity)
    # Add similarity scores as a new column in the DataFrame
    df['similarity'] = similarity_scores
    # Sort DataFrame by sort_column and similarity
    df_sorted = df.sort_values(by='similarity', ascending=False)
    # Return the IDs from each row as a list
    return df_sorted['id'].tolist()

