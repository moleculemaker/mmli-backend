import json

from services.minio_service import MinIOService

class OEDService:
    
    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService):
        file = service.get_file(bucket_name, f"{job_id}/out/output.json")
        if not file:
            return None
        data = json.loads(file)

        return data