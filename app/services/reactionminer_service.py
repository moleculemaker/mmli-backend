import json
import os
import time

from fastapi import HTTPException
from sqlmodel.ext.asyncio.session import AsyncSession

from config import log
from models.enums import JobStatus
from services.minio_service import MinIOService

from services.reactionminer_search.core import txt_eval, smi_eval, mm_eval


class ReactionMinerService:
    def __init__(self, db) -> None:
        self.db = db

    async def update_job_phase(self, jobObject, phase: JobStatus):
        jobObject.phase = phase
        if phase == JobStatus.PROCESSING:
            jobObject.time_start = int(time.time())
        else:
            jobObject.time_end = int(time.time())
        self.db.add(jobObject)
        await self.db.commit()

    @staticmethod
    async def resultPostProcess(bucket_name: str, job_id: str, service: MinIOService, db: AsyncSession):
        """
        Inputs stored in Minio:  /{job_id}/in/[name].pdf    Bucket name: reactionminer
        Outputs stored in Minio: /{job_id}/out/[name].json  Bucket name: reactionminer
        """
        folder_path = f"/{job_id}/out/"
        objects = service.list_files(bucket_name, folder_path, True)

        # Iterate over folder and add all contents to a dictionary
        content = {}
        for obj in objects:
            file_name = os.path.basename(obj.object_name).split('/')[-1]
            if file_name.endswith('.json'):
                content[file_name] = json.loads(service.get_file(bucket_name=bucket_name, object_name=obj.object_name))
            elif file_name.endswith('.csv'):
                content[file_name] = service.get_file(bucket_name=bucket_name, object_name=obj.object_name)
            else:
                log.warning(f'Skipping unrecognized file extension: ' + str(file_name))

        # Return the dictionary if it has contents
        if not content:
            raise HTTPException(status_code=404, detail=f"No output files were found")

        return content

    @staticmethod
    # XXX: we need to keep this method to maintain the reactionminer-demo instance for SIGIR
    # See https://reactionminer-demo.platform.moleculemaker.org/reaction-miner
    def search_reactionminer(smi_search: str, text_search: str):
        root_path = os.path.dirname(os.path.abspath(__file__))
        par_path = os.path.join(root_path, 'reactionminer_search/indexing_data_suzuki/tsv_files_full/Pars.tsv')
        smi_path = os.path.join(root_path, 'reactionminer_search/indexing_data_suzuki/tsv_files_full/Par_SMI.tsv')
        par2smi_path = os.path.join(root_path, 'reactionminer_search/indexing_data_suzuki/tsv_files_full/Pars2SMI.tsv')
        par2img_path = os.path.join(root_path, 'reactionminer_search/indexing_data_suzuki/tsv_files_full/Par_Imgs.tsv')
        smi2name_path = os.path.join(root_path, 'reactionminer_search/indexing_data_suzuki/tsv_files_full/Par_SMI_Names.tsv')

        try:
            # Text Only Search
            if not len(smi_search):
                result, img_result = txt_eval(par_path, par2smi_path, par2img_path,
                                smi2name_path, smi_path, txt_q=text_search, topk=10)
            # SMILES Only Search
            elif not len(text_search):
                result, img_result = smi_eval(par_path, smi_path, par2smi_path,
                                par2img_path, smi2name_path,
                                choice='tani', smi_q=smi_search)
            # Multi-modal Search
            else:
                result, img_result = mm_eval(par_path, smi_path, par2smi_path,
                                par2img_path, smi2name_path, txt_q=text_search, smi_q=smi_search,
                                choice='sub', topk=10, fm='priority')

            result_json = result.to_dict(orient='records')
            img_result_json = img_result.to_dict(orient='records')

            return {
                'status': 'success',
                'passage_results': result_json,
                'image_results': img_result_json
            }

        # the code throws a KeyError if the search SMILES is invalid
        # in this case, we return an empty list
        except KeyError as e:
            return {
                'status': 'success',
                'passage_results': [],
                'image_results': []
            }

        # for all other exceptions, we return the error
        except Exception as e:
            return {
                'status': 'error',
                'error': str(e)
            }