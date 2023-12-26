from typing import List
from pydantic import BaseModel

class AnalyzeRequestBody(BaseModel):
    jobId: str
    user_email: str
    fileList: List[str]
