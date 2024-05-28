import json

from config import get_logger

log = get_logger(__name__)

class MolliService:

    @staticmethod
    def build_molli_job_environment(job_info, job_id):
        environment = []

        input_dir = f'/uws/jobs/molli/{job_id}/in'
        CORES_FILE_NAME = job_info['CORES_FILE_NAME']
        SUBS_FILE_NAME = job_info['SUBS_FILE_NAME']

        # Tell the job to use these files when it runs
        environment.append({'name': 'CORES_INPUT_FILE', 'value': f'{input_dir}/{CORES_FILE_NAME}'})
        environment.append({'name': 'SUBS_INPUT_FILE', 'value': f'{input_dir}/{SUBS_FILE_NAME}'})

        return environment

    @staticmethod
    async def molliResultPostProcess(bucket_name, job_id, service, db):
        gen = service.get_file(bucket_name, f"{job_id}/out/test_combine_new_env_library.json")
        pca = service.get_file(bucket_name, f"{job_id}/out/new_env_data3_pca")
        tsne = service.get_file(bucket_name, f"{job_id}/out/new_env_data3_tsne")

        result = {
            'structures': json.loads(gen),
            'clusteringData': {
                'pca': json.loads(pca),
                'tsne': json.loads(tsne)
            }
        }

        return json.dumps(result)
