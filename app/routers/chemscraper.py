
import csv
import uuid
import io
from typing import List
from datetime import datetime
from fastapi import APIRouter, Depends, HTTPException, UploadFile, File, status, Query, BackgroundTasks
from fastapi.responses import JSONResponse, FileResponse
import zipfile

from services.minio_service import MinIOService

from services.chemscraper_service import ChemScraperService

from models.analyzeRequestBody import AnalyzeRequestBody

from models.exportRequestBody import ExportRequestBody
from typing import Optional


import pandas as pd

router = APIRouter()


@router.post("/chemscraper/analyze", tags=['ChemScraper'])
async def analyze_documents(requestBody: AnalyzeRequestBody, background_tasks: BackgroundTasks, service: MinIOService = Depends()):
    # Analyze only one document for NSF demo
    if len(requestBody.fileList) > 0 and requestBody.jobId != "":
        filename = requestBody.fileList[0]
        chemscraperService = ChemScraperService()
        objectPath = f"inputs/{requestBody.jobId}/{filename}"
        background_tasks.add_task(chemscraperService.runChemscraperOnDocument, 'chemscraper', filename, objectPath, requestBody.jobId, service)
        content = {"jobId": requestBody.jobId, "submitted_at": datetime.now().isoformat()}
        return JSONResponse(content=content, status_code=status.HTTP_202_ACCEPTED) 
