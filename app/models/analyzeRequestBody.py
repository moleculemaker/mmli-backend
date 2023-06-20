from typing import List
from pydantic import BaseModel

class AnalyzeRequestBody(BaseModel):
    jobId: str
    user_email: str
    captcha_token: str
    fileList: List[str]
