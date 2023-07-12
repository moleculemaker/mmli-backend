from fastapi import HTTPException
import requests
import io
import os
from services.minio_service import MinIOService
from fastapi.responses import JSONResponse

class ChemScraperService:
    chemscraper_api_baseURL = os.environ.get("CHEMSCRAPER_API_BASE_URL")

    def __init__(self) -> None:
        pass

    def runChemscraperOnDocument(self, bucket_name: str, filename: str, objectPath: str, jobId: str, service: MinIOService):
        data = service.get_file(bucket_name, objectPath)
        if data is None:
            raise HTTPException(status_code=404, detail="File not found")
        
        data_bytes = io.BytesIO(data)
        response = requests.post(self.chemscraper_api_baseURL + '/extractPdf', files={'pdf': (filename, data_bytes)})

        if response.status_code == 200:
            tsv_content = response.text.encode()
            upload_result = service.upload_file(bucket_name, "results/" + jobId + ".tsv", tsv_content)
            if upload_result:
                return True
        else:
            error_content = response.text.encode()
            upload_result = service.upload_file(bucket_name, "errors/" + jobId + ".txt", error_content)
        return False