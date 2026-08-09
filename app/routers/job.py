import base64
import os
import json
import re
import time
import uuid
import csv
import io
import traceback

from typing import List

from fastapi import Depends, HTTPException, APIRouter, UploadFile
from fastapi.openapi.models import Response
from fastapi.params import Path, Body, File
from pydantic.fields import Annotated, Optional
from sqlalchemy import delete
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession
from starlette import status
from starlette.responses import JSONResponse

from config import get_logger, app_config, STATUS_ERROR
from models.enums import JobType, JobStatus, JobTypes
from models.sqlmodel.db import get_session
from models.sqlmodel.models import Job, JobCreate, JobUpdate

from services import kubejob_service
from services.clean_service import CleanService
from services.crispr_copies_service import CRISPRCopiesService
from services.molli_service import MolliService
from services.minio_service import MinIOService
from services.mutagenesis_service import MutagenesisService
from services.shared import is_valid_pdb_file
from services.somn_service import SomnService

router = APIRouter()

log = get_logger(__name__)

# EZSpecificity input limits (mirror the frontend's MAX_ENZYMES / substrate cap)
EZSPEC_MAX_ENZYMES = 5
EZSPEC_MAX_SUBSTRATES = 10


async def _mark_job_failed(db: AsyncSession, db_job: Job) -> None:
    """Move a just-created job to 'error' when nothing will ever run to advance it.

    Only ever called for a row this request created. A job that already existed must
    keep whatever phase it legitimately reached.
    """
    db_job.phase = JobStatus.ERROR
    db.add(db_job)
    await db.commit()


CREATE_JOB_DESCRIPTION = """
Create a new run for a new or existing Job.

`job_type` is the path segment, e.g. `crispr-copies`, `mutagenesis`, `clean`, `somn`.

`job_info` is a **JSON-encoded string** (not a nested object; future versions may change this)
whose schema depends on `job_type`. For jobs that take uploaded files, first upload each
file via `POST /{job_type}/upload?job_id=<id>` (using the same `job_id`), then reference
the uploaded filenames here.

Results are retrieved from `GET /{job_type}/results/{job_id}` once the job completes.
"""
@router.post("/{job_type}/jobs", response_model=Job, tags=['Jobs'],
             summary="Create a job run", description=CREATE_JOB_DESCRIPTION)
async def create_job(
        job_id: Optional[str] = Body(default=None),
        run_id: Optional[str] = Body(default=None),
        email: Optional[str] = Body(default=None),
        job_info: Optional[str] = Body(default="{}"),
        job_type: str = Path(),
        service: MinIOService = Depends(),
        db: AsyncSession = Depends(get_session)
):
    job_id = job_id if job_id else str(kubejob_service.generate_uuid())
    #run_id = run_id if run_id else str(kubejob_service.generate_uuid())

    # Check if this job_id already exists for this job_type
    statement = select(Job).where(Job.type == job_type).where(Job.job_id == job_id)
    existing_jobs = await db.exec(statement)
    db_job: Job = existing_jobs.first()
    #if db_job:
    #    raise HTTPException(status_code=409, detail=f"Job already exists with job_id={job_id}")

    # Validate Job type
    # TODO: Set command+image based on job_type
    if job_type in JobTypes:
        log.debug(f"Creating Kubernetes job: {job_type}")
        # runs in Kubernetes, read Docker image name from config
        image_name = app_config['kubernetes_jobs'][job_type]['image']

        # Command + environment are set differently for each job (see below)
        command = ''
        environment = []

        # Mount in secrets/volumes at runtime
        #volumes = []
        #secrets = []

        if job_type == JobType.DEFAULT:
            command = app_config['kubernetes_jobs'][job_type]['command']
            environment = app_config['kubernetes_jobs'][job_type]['env'] if 'env' in app_config['kubernetes_jobs'][job_type] else []
            #initContainers = app_config['kubernetes_jobs'][job_type]['initContainers'] if 'initContainers' in app_config['kubernetes_jobs'][job_type] else []
            #image = app_config['kubernetes_jobs'][job_type]['image'] if 'image' in app_config['kubernetes_jobs'][job_type] else app_config['kubernetes_jobs']['defaults']
            #imagePullPolicy = app_config['kubernetes_jobs'][job_type]['imagePullPolicy'] if 'imagePullPolicy' in app_config['kubernetes_jobs'][job_type] else 'Always'
            #volumes = app_config['kubernetes_jobs'][job_type]['volumes'] if 'volumes' in app_config['kubernetes_jobs'][job_type] else []
            #secrets = app_config['kubernetes_jobs'][job_type]['secrets'] if 'secrets' in app_config['kubernetes_jobs'][job_type] else []

        elif job_type == JobType.ACERETRO:
            # ACERetro jobs
            # Example usage:
            # curl -X POST https://mmli.kastan.ai/aceretro/jobs \
            #   -H "Content-Type: application/json" \
            #   -d '{
            #     "job_id": "123",
            #     "email": "user@gmail.com",
            #     "job_info": "{\"smiles\": \"O=C(COP(=O)(O)O)[C@H](O)[C@H](O)CO\"}"
            #   }'
            log.info(f"------------------ STARTING ACERETRO JOB ------------------  job[{job_type}]: " + job_id)
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            
            from config import MINIO_SERVER, MINIO_ACCESS_KEY, MINIO_SECRET_KEY
            environment = [{
                'name': 'MINIO_URL',
                # 'value': app_config['minio']['apiBaseUrl']
                'value': MINIO_SERVER
            },
            {
                'name': 'MINIO_ACCESS_KEY',
                # 'value': app_config['minio']['accessKey']
                'value': MINIO_ACCESS_KEY
            },
            {
                'name': 'MINIO_SECRET_ACCESS_KEY',
                # 'value': app_config['minio']['accessKey']
                'value': MINIO_SECRET_KEY
            },
            {
                'name': 'MINIO_SECURE',
                # 'value': app_config['minio']['secretKey']
                'value': False # HARD CODED, like everywhere else in this repo
            }]

            command = f"python entrypoint.py --job_id {job_id}"
            # Job is created at end of function

        elif job_type == JobType.CHEMSCRAPER:
            log.debug(f"Creating Kubernetes job: {job_type}")

            job_config = json.loads(job_info.replace('\"', '"'))
            if 'input_file' not in job_config:
                raise HTTPException(status_code=400, detail='"job_info" requires "input_file" for ChemScraper jobs')

            environment = [
                {
                    'name': 'CHEMSCRAPER_INPUT_FILE',
                    'value': job_config['input_file']
                }
            ]

        elif job_type == JobType.CLEAN:
            # Build up input.FASTA from user input
            job_config = json.loads(job_info.replace('\"', '"'))
            command = CleanService.build_clean_job_command(job_id=job_id, job_info=job_config)

        elif job_type == JobType.CLEANDB_MEPESM:
            # Example: "job_info": "{\"sequence\":\"MEDIPDTSRPPLKYVK...\"}"
            # OR       "job_info": "MEDIPDTSRPPLKYVK..."
            try:
                log.debug(f'    job_info: {job_info}')
                job_config = json.loads(job_info.replace('\\"', '"'))
                log.debug(f'    job_config: {job_config}')
                input_sequence = job_config['sequence']
                log.debug(f'CLEANDB-MEPESM - Parsed job_info as JSON')
                environment = [{'name': 'CLEANDB_INPUT_SEQUENCE', 'value': input_sequence}]
            except Exception as e:
                log.warning(f'WARNING: CLEANDB-MEPESM - Failed to parse job_info as JSON: {str(e)}')
                traceback.format_exc()
                log.warning(f'WARNING: CLEANDB-MEPESM - Falling back to treating job_info as a string')
                environment = [{'name': 'CLEANDB_INPUT_SEQUENCE', 'value': job_info}]

            log.debug(f'    environment: {environment}')

        elif job_type == JobType.CRISPR_COPIES:
            # Inputs come from either uploaded files (synced from MinIO {job_id}/in/ by
            # prejob.py) or a dropdown organism accession that the init container's
            # fetch_organism.py resolves into ${JOB_INPUT_DIR}. Pass the accession (if the
            # dropdown path was used) so the init container knows what to fetch.
            job_config = json.loads(job_info.replace('\"', '"'))
            organism_accession = (
                (job_config.get("organism") or {}).get("organismIdentifier") or {}
            ).get("value") or ""
            environment = [{'name': 'ORGANISM_ACCESSION', 'value': organism_accession}]
            command = CRISPRCopiesService.build_crispr_copies_job_command(job_id=job_id, job_info=job_config)

        # EZspecificity parent job, see subjobs below
        elif job_type == JobType.EZ_SPECIFICITY:
            command = app_config['kubernetes_jobs'][job_type]['command']
            job_config = json.loads(job_info.replace('\"', '"'))

            # Validate user input: we need enzymes + substrates to build the docking config
            if 'enzymes' not in job_config or 'substrates' not in job_config:
                raise HTTPException(status_code=400,
                    detail='"job_info" requires "enzymes" and "substrates"')

            enzymes = job_config['enzymes']
            substrates = job_config['substrates']
            if not isinstance(enzymes, list) or not isinstance(substrates, list):
                raise HTTPException(status_code=400,
                    detail='"enzymes" and "substrates" must be lists')
            if not 1 <= len(enzymes) <= EZSPEC_MAX_ENZYMES:
                raise HTTPException(status_code=400,
                    detail=f'Expected 1-{EZSPEC_MAX_ENZYMES} enzyme(s), got {len(enzymes)}')
            if not 1 <= len(substrates) <= EZSPEC_MAX_SUBSTRATES:
                raise HTTPException(status_code=400,
                    detail=f'Expected 1-{EZSPEC_MAX_SUBSTRATES} substrate(s), got {len(substrates)}')

            # Validate each uploaded enzyme structure. The PDBs were already uploaded
            # to {job_id}/in/ via the /upload endpoint; we re-check them here (server-side,
            # where we can actually parse the file content) before kicking off the pipeline.
            for enzyme in enzymes:
                filename = enzyme.get('filename')
                if not filename:
                    raise HTTPException(status_code=400, detail='Each enzyme requires a "filename"')
                content = service.get_file(job_type, f"{job_id}/in/{filename}")
                if content is None:
                    raise HTTPException(status_code=400,
                        detail=f'Enzyme structure not found in uploads: {filename}')
                if not is_valid_pdb_file(content):
                    raise HTTPException(status_code=400,
                        detail=f'Invalid PDB structure: {filename}')

            # Write job_config.json into the parent job's input dir. coordinator.py copies
            # the parent's in/ into each subjob, so unidock + inference both receive this
            # config (the containers read job_config.json, not the job_info API field).
            container_config = {
                'enzymes': enzymes,
                'substrates': substrates,
            }
            if 'docking' in job_config:
                container_config['docking'] = job_config['docking']
            if service.ensure_bucket_exists(job_type):
                service.upload_file(
                    job_type,
                    f"{job_id}/in/job_config.json",
                    json.dumps(container_config).encode('utf-8'),
                )

            # Generate subjob_ids, store these in job_info / pass along as environment
            ezspec_unidock_job_id = str(kubejob_service.generate_uuid())
            ezspec_inference_job_id = str(kubejob_service.generate_uuid())

            # Preserve these subjob_ids in our job_info
            job_config = job_config | {
                'parent_job_id': job_id,
                'ezspec_unidock_job_id': ezspec_unidock_job_id,
                'ezspec_inference_job_id': ezspec_inference_job_id,
            }
            job_info = json.dumps(job_config).replace('"', '\"')

            # Pass parent / subjob_ids along as envvars
            environment += app_config['kubernetes_jobs'][job_type]['env']
            environment += [
                { "name": "PARENT_JOB_ID", "value": job_id },
                { "name": "EZSPEC_UNIDOCK_JOB_ID", "value": ezspec_unidock_job_id },
                { "name": "EZSPEC_INFERENCE_JOB_ID", "value": ezspec_inference_job_id }
            ]

        # all EZspecificity subjobs / job steps share the same handling
        elif job_type == JobType.EZSPEC_UNIDOCK or job_type == JobType.EZSPEC_INFERENCE:
            # TODO: update command to handle ez-specificity jobs
            #command = app_config['kubernetes_jobs'][job_type]['command']

            # Grab our subjob_ids from the passed job_info
            job_config = json.loads(job_info.replace('\"', '"'))
            if job_type == JobType.EZSPEC_UNIDOCK:
                job_id = job_config['ezspec_unidock_job_id']
            elif job_type == JobType.EZSPEC_INFERENCE:
                job_id = job_config['ezspec_inference_job_id']

            # Pass parent / subjob_ids along as envvars
            environment += app_config['kubernetes_jobs'][job_type]['env']
            environment += [
                { "name": "PARENT_JOB_ID",  "value": job_config['parent_job_id'] },
                { "name": "EZSPEC_UNIDOCK_JOB_ID", "value": job_config['ezspec_unidock_job_id'] },
                { "name": "EZSPEC_INFERENCE_JOB_ID", "value": job_config['ezspec_inference_job_id'] }
            ]

        elif job_type == JobType.ML_SIMPLEFOLD:
            log.info(f"------------------ STARTING ML-SIMPLEFOLD JOB ------------------  job[{job_type}]: " + job_id)
            job_config = json.loads(job_info.replace('\"', '"'))

            if 'fasta' not in job_config:
                raise HTTPException(status_code=400, detail='"job_info" requires "fasta" for SimpleFold jobs')

            # Upload FASTA content to MinIO
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.fasta", job_config['fasta'].encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload FASTA to MinIO")

            command = (
                "simplefold"
                " --simplefold_model simplefold_100M"
                " --num_steps 500"
                " --tau 0.01"
                " --nsample_per_protein 1"
                " --plddt"
                " --fasta_path ${JOB_INPUT_DIR}/input.fasta"
                " --output_dir ${JOB_OUTPUT_DIR}"
                " --backend torch"
                " && rm -rf ${JOB_OUTPUT_DIR}/cache"
            )

        elif job_type == JobType.MOLLI:
            # Pass path to CORES/SUBS files into the container
            command = app_config['kubernetes_jobs'][job_type]['command']
            job_config = json.loads(job_info.replace('\"', '"'))
            environment = MolliService.build_molli_job_environment(job_id=job_id, job_info=job_config)

        elif job_type == JobType.MUTAGENESIS:
            # Inputs (mutation list / ORF file) are uploaded by the frontend to
            # MinIO {job_id}/in/ and synced into ${JOB_INPUT_DIR} before the job runs.
            job_config = json.loads(job_info.replace('\"', '"'))
            command = MutagenesisService.build_mutagenesis_job_command(job_id=job_id, job_info=job_config)

        elif job_type == JobType.NOVOSTOIC_DGPREDICTOR:
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']

        elif job_type == JobType.NOVOSTOIC_ENZRANK:
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']

        elif job_type == JobType.NOVOSTOIC_OPTSTOIC:
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']

            # environment = [{
            #     # TBD... 
            #     # 'name': 'SOMN_PROJECT_DIR',
            #     # 'value': somn_project_dir
            # }]

        elif job_type == JobType.NOVOSTOIC_PATHWAYS:
            if service.ensure_bucket_exists(job_type):
                job_info = json.loads(job_info.replace('\"', '"'))
                stoic = f'{job_info["substrate"]["amount"]} {job_info["substrate"]["molecule"]}'
                for coReactant in job_info['reactants']:
                    stoic += f' + {coReactant["amount"]} {coReactant["molecule"]}'
                stoic += " <=> "
                for coProduct in job_info['products']:
                    stoic += f'{coProduct["amount"]} {coProduct["molecule"]} + '
                stoic += f'{job_info["product"]["amount"]} {job_info["product"]["molecule"]}'
                
                job_info['stoic'] = stoic
                job_info['substrate'] = job_info['substrate']['molecule']
                job_info['product'] = job_info['product']['molecule']
                job_info['num_enzymes'] = job_info['num_enzymes'] if 'num_enzymes' in job_info else 0
                
                job_info = json.dumps(job_info)
                upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']

        elif job_type == JobType.OED_CHEMINFO:
            # Pass path to CORES/SUBS files into the container
            if service.ensure_bucket_exists(job_type):
                upload_result = service.upload_file(job_type, f"/{job_id}/in/job.json", job_info.replace('\"', '"').encode('utf-8'))
                if not upload_result:
                    raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
            command = app_config['kubernetes_jobs'][job_type]['command']

        elif job_type == JobType.OED_DLKCAT or job_type == JobType.OED_UNIKP or job_type == JobType.OED_CATPRED:
            # Example: "job_info": "{\"input_pairs\":[{\"name\":\"example\",\"sequence\":\"MEDIPDTSRPPLKYVK...\",\"type\":\"FASTA\",\"smiles\":\"OC1=CC=C(C[C@@H](C(O)=O)N)C=C1\"}]}"
            log.debug(f'Running OpenEnzymeDB job: {job_type} - {job_id}')
            log.debug(f'    job_info: {job_info}')
            job_config = json.loads(job_info.replace('\\"', '"'))
            log.debug(f'    job_config: {job_config}')
            job_config_str = json.dumps(job_config['input_pairs']).replace('"', '\\"')
            environment = [{'name': 'OED_INPUT_PAIRS', 'value': job_info.replace('"', '\\"')}]
            log.debug(f'    environment: {environment}')

        elif job_type == JobType.REACTIONMINER:
            log.debug(f'Running ReactionMiner job: {job_id}')
            environment = app_config['kubernetes_jobs']['reactionminer']['env']

        elif job_type == JobType.SOMN:
            #  Build up example_request.csv from user input, upload to MinIO?
            json_str = job_info.replace('\"', '"')
            job_config = json.loads(json_str)
            
            # Canonicalize SMILES and update names from reference files
            for config in job_config:
                if config['el_input_type'] == 'smi':
                    config['el'] = SomnService.canonicalize_smiles(config['el'])
                if config['nuc_input_type'] == 'smi':
                    config['nuc'] = SomnService.canonicalize_smiles(config['nuc'])
            
            # Generate unique name mappings
            job_config, el_name_map = SomnService.generate_name_mapping(job_config, 'el')
            job_config, nuc_name_map = SomnService.generate_name_mapping(job_config, 'nuc')
            
            # job_config, el_name_map = SomnService.update_names_from_reference(job_config, el_name_map, 'el')
            # job_config, nuc_name_map = SomnService.update_names_from_reference(job_config, nuc_name_map, 'nuc')
            
            print('updated config: ', job_config, el_name_map, nuc_name_map)
            
            # job_info is automatically stored in Postgres to retain user input
            job_info = json.dumps({
                'info': job_config,
                'el_name_map': el_name_map,
                'nuc_name_map': nuc_name_map
            })
            
            # Some SMILESs fail to pass 3d generation test, so we process them as mol2 input
            for config in job_config:
                config['el'] = SomnService.process_molecule_input(config, 'el')
                config['nuc'] = SomnService.process_molecule_input(config, 'nuc')
            
            
            file = io.StringIO()
            writer = csv.writer(file)
            writer.writerow([
                "user", 
                "nuc", 
                "el", 
                "nuc_name", 
                "el_name",
                "nuc_idx",
                "el_idx"
            ])
            for config in job_config:
                writer.writerow([
                    config["reactant_pair_name"], 
                    config["nuc"],
                    config["el"], 
                    config["nuc_name"], 
                    config["el_name"],
                    config["nuc_idx"],
                    config["el_idx"]
                ])
            
            upload_result = service.upload_file(job_type, f"/{job_id}/in/example_request.csv", file.getvalue().encode('utf-8'))
            if not upload_result:
                raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")

            # We assume that file has already been uploaded to MinIO
            somn_project_dir = app_config['kubernetes_jobs']['somn']['projectDirectory']

            command = app_config['kubernetes_jobs'][job_type]['command']
            #command = f"ls -al && whoami && somn predict {project_id} {model_set} {new_predictions_name}"

            environment = [{
                'name': 'SOMN_PROJECT_DIR',
                'value': somn_project_dir
            }]

        # TODO: Set user_agent based on requestor
        user_agent = ''

        # Write the DB row BEFORE creating the Kubernetes Job. KubeWatcher looks this row
        # up by jobId for every event it receives and drops events whose row does not
        # exist yet; since the watch never replays old events, a Job that got created
        # (and, for the short ones, finished) inside that window is stranded at 'queued'
        # even after Kubernetes reports it complete.
        #
        # The same reordering exists on the open PR that fixes that stranding bug. It is
        # repeated here rather than depended upon, because the duplicate-submission
        # handling below needs the row to already exist and this branch does not build on
        # that PR. Expect a small conflict if both land.
        created_new_job = not db_job
        if created_new_job:
            db_job: Job = Job(
                # User input
                email=email,
                job_info=job_info,
                job_id=job_id,
                run_id=run_id,

                type=job_type,
                command=command,
                image=image_name,

                # Job metadata
                deleted=0,
                time_created=int(time.time()),

                # Set ser metadata
                user_agent=user_agent,
            )

            db.add(db_job)
            await db.commit()
            await db.refresh(db_job)

        # A POST naming a job_id that already exists is not a resubmission, and this
        # endpoint has no authentication: the caller has proved only that they can name
        # an id, which is not authority over the job that holds it. Stop here.
        #
        # Previously this path fell through and (a) launched a second Kubernetes Job
        # under the same name, (b) replaced the stored job_info with the caller's, and
        # (c) returned the owner's email address and inputs in the response body.
        # The response keys and status code are unchanged so existing clients - notably
        # coordinator.py, which reads job_id from a 200 or a 201 - keep working.
        if not created_new_job:
            log.warning(f'Ignoring duplicate job submission for existing job_id={job_id} (type={job_type})')
            return JSONResponse(status_code=status.HTTP_200_OK, content={
                'job_id': str(db_job.job_id),
                'run_id': str(db_job.run_id),
                # Withheld: naming an id is not authorisation to read its owner's
                # address or the inputs they submitted.
                'email': None,
                'job_info': None,
            })

        # Run a Kubernetes Job with the given image + command + environment
        try:
            log.debug(f"Creating Kubernetes job[{job_type}]: " + job_id)
            create_response = kubejob_service.create_job(
                job_type=job_type,
                job_id=job_id,
                run_id=run_id,
                image_name=image_name,
                command=command,
                environment=environment
            )
        except Exception as ex:
            log.error("Failed to create Job: " + str(ex))
            log.error(traceback.format_exc())
            # The row we just wrote is already visible to clients, and nothing will ever run
            # to move it along - fail it here rather than leaving it polling 'queued' forever.
            await _mark_job_failed(db, db_job)
            raise HTTPException(status_code=400, detail="Failed to create Job: " + str(ex))

        # create_job catches ApiException internally, logs it, and returns a response
        # carrying status=error rather than propagating it. Kubernetes rejections -
        # an exceeded quota, an unpullable image, a name that already exists - are
        # therefore invisible to the except clause above. Without this check the caller
        # is told 201 while the row sits at 'queued' forever with no pod behind it.
        if isinstance(create_response, dict) and create_response.get('status') == STATUS_ERROR:
            message = create_response.get('message', 'unknown error')
            log.error(f'Kubernetes rejected job[{job_type}] {job_id}: {message}')
            await _mark_job_failed(db, db_job)
            raise HTTPException(status_code=400, detail="Failed to create Job: " + str(message))

        return JSONResponse(status_code=status.HTTP_201_CREATED, content={
            'job_id': str(db_job.job_id),
            'run_id': str(db_job.run_id),
            'email': str(db_job.email),
            'job_info': str(db_job.job_info),
        })

    else:
        log.error("Failed to create job - invalid job type: " + str(job_type))
        raise HTTPException(status_code=400, detail="Invalid job type: " + str(job_type))


@router.get("/{job_type}/jobs", tags=['Jobs'], description="Get a list of all job runs by type")
async def list_jobs_by_type(job_type: str, db: AsyncSession = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    results = await db.execute(select(Job).where(Job.type == job_type))
    jobs = results.scalars().all()
    return jobs


@router.get("/{job_type}/jobs/{job_id}", tags=['Jobs'], description="Get a list of all job runs by type and job_id")
async def list_jobs_by_type_and_job_id(job_type: str, job_id: str, db: AsyncSession = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    results = await db.execute(select(Job).where(Job.type == job_type).where(Job.job_id == job_id))
    jobs = results.scalars().all()
    return jobs


@router.get("/{job_type}/jobs/{job_id}/{run_id}", tags=['Jobs'], description="Get a single job by type, job_id, and run_id")
async def get_job_by_type_and_job_id_and_run_id(job_type: str, job_id: str, run_id: str, db: AsyncSession = Depends(get_session)):
    job_types = [e for e in JobType]
    if job_type not in job_types:
        raise HTTPException(status_code=400, detail="Invalid job type: " + job_type)

    result = await db.execute(select(Job).where(Job.job_id == job_id).where(Job.run_id == run_id))
    job = result.first()
    if not job:
        raise HTTPException(status_code=404, detail=f"No job found with job_id={job_id} and run_id={run_id}")
    return job.Job


@router.put("/{job_type}/jobs/{job_id}/{run_id}", response_model=Job, tags=['Jobs'], description="Overwrite all writeable fields of an existing Job")
async def update_existing_job(job: Job, job_type: str, db: AsyncSession = Depends(get_session)):
    # Check if this job_id already exists
    db_job = await db.get(Job, job.job_id)
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job.job_id} and run_id={job.run_id}")

    # TODO: Validation
    # TODO: Set internal metadata / fields
    db_job.phase = job.phase
    db_job.time_start = job.time_start
    db_job.time_end = job.time_end

    # Update existing DB job from user input
    db.add(db_job)
    await db.commit()
    await db.refresh(db_job)

    return db_job


@router.patch("/{job_type}/jobs/{job_id}/{run_id}", response_model=Job, tags=['Jobs'], description="Update one or more fields of an existing Job")
async def patch_existing_job(job: JobUpdate, job_type: str, db: AsyncSession = Depends(get_session)):
    # Check if this job_id already exists
    db_job = await db.get(Job, job.job_id)
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job.job_id} and run_id={job.run_id}")

    # TODO: Validation
    # TODO: Set internal metadata / fields
    db_job.phase = job.phase if job.phase is not None else db_job.phase
    db_job.time_start = job.time_start if job.time_start is not None else db_job.time_start
    db_job.time_end = job.time_end if job.time_end is not None else db_job.time_end

    # Update existing DB job from user input
    db.add(db_job)
    await db.commit()
    await db.refresh(db_job)

    return db_job


@router.delete("/{job_type}/jobs/{job_id}/{run_id}", tags=['Jobs'], description="Delete a single Job by type, job_id, and run_id")
async def delete_job_by_type_and_job_id_and_run_id(job_type: str, job_id: str, run_id: str, db: AsyncSession = Depends(get_session)):
    # Check if this job_id already exists
    result = await db.execute(select(Job).where(Job.type == job_type).where(Job.job_id == job_id).where(Job.run_id == run_id))
    db_job = result.first()
    if not db_job:
        raise HTTPException(status_code=404, detail=f"Job does not exist with type={job_type} job_id={job_id} and run_id={run_id}")

    # Delete the Kubernetes Job before dropping the row. Removing only the row leaves a
    # pod running with nothing tracking it: KubeWatcher looks each event up by jobId and
    # discards events whose row is gone, so the work continues, consumes cluster
    # resources, and writes its output to MinIO with no way to observe or stop it.
    #
    # Deliberately best-effort. delete_job already swallows ApiException (a Job that
    # Kubernetes has since reaped is a 404, which is the expected case for anything
    # older than ttlSecondsAfterFinished), and a cluster problem must not leave the
    # caller unable to delete their own record.
    try:
        kubejob_service.delete_job(job_type=job_type, job_id=job_id)
    except Exception as ex:
        log.error(f'Failed to delete Kubernetes job[{job_type}] {job_id}, deleting DB row anyway: {ex}')

    await db.execute(delete(Job).where(Job.type == job_type).where(Job.job_id == job_id).where(Job.run_id == run_id))
    await db.commit()

    return db_job.Job
