from config import get_logger

log = get_logger(__name__)

class EzSpecificityService:

    @staticmethod
    async def resultPostProcess(bucket_name, job_id, service, db):
        # TODO: what format is expected for the results in the frontend?

        # TODO: fetch MinIO files and present in the correct format
        enzymes = service.get_file(bucket_name, f"{job_id}/out/Enzymes.csv")
        substrates = service.get_file(bucket_name, f"{job_id}/out/Substrates.csv")
        #data = service.get_file(bucket_name, f"{job_id}/out/data.csv")
        #job_config = service.get_file(bucket_name, f"{job_id}/out/job_config.json")

        return {
            'hello': 'world'
        }

    @staticmethod
    async def unidockResultPostProcess(bucket_name, job_id, service, db):
        # TODO: what format is expected for the results in the frontend?

        # TODO: fetch MinIO files and present in the correct format

        return {
            'hello': 'unidock'
        }

    @staticmethod
    async def inferenceResultPostProcess(bucket_name, job_id, service, db):
        # TODO: what format is expected for the results in the frontend?

        # TODO: fetch MinIO files and present in the correct format

        return {
            'hello': 'inference'
        }
