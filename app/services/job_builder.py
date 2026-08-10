"""Per-tool job preparation.

Turns a job type and its `job_info` into the image, command and environment the
Kubernetes job needs, plus any side effects a tool requires before it can run -- most
often writing an input file to MinIO.

This lives outside the router because two APIs need it: the legacy `/{job_type}/jobs`
endpoint and the versioned `/v1` one. It is a straight extraction, deliberately not a
rewrite: the dispatch below is the same chain that has always run, moved rather than
reworked, so that the characterization tests keep their meaning. Tidying it is a
separate job from making it reusable.

Two tools rewrite what they were given, which is why both `job_id` and `job_info` come
back out:

  novostoic-pathways  assembles the `stoic` string its container expects and flattens
                      substrate/product, storing the derived form
  somn                canonicalizes SMILES and de-duplicates names, storing the mapping
                      alongside the user's original input
  ez-specificity      mints subjob ids and records them in job_info
  ezspec-*            take their own job_id from the parent's job_info
"""
import csv
import io
import json
import traceback
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from fastapi import HTTPException

from config import get_logger, app_config
from models.enums import JobType
from services import kubejob_service
from services.clean_service import CleanService
from services.minio_service import MinIOService
from services.molli_service import MolliService
from services.shared import is_valid_pdb_file
from services.somn_service import SomnService

log = get_logger(__name__)

# EZSpecificity input limits (mirror the frontend's MAX_ENZYMES / substrate cap)
EZSPEC_MAX_ENZYMES = 5
EZSPEC_MAX_SUBSTRATES = 10


@dataclass
class PreparedJob:
    """Everything needed to launch a job, after per-tool preparation."""
    job_id: str
    job_info: str
    image_name: str
    command: str = ''
    environment: List[Dict[str, Any]] = field(default_factory=list)
    parent_job_id: Optional[str] = None


def prepare_job(job_type: str, job_id: str, job_info: str, service: MinIOService) -> PreparedJob:
    """Run a tool's pre-launch preparation and return what the job needs to run.

    Raises HTTPException for invalid input, exactly as the inline version did.
    """
    log.debug(f"Creating Kubernetes job: {job_type}")
    # runs in Kubernetes, read Docker image name from config
    image_name = app_config['kubernetes_jobs'][job_type]['image']

    # Command + environment are set differently for each job (see below)
    command = ''
    environment = []
    # Only set for subjobs created by a parent job's coordinator
    parent_job_id = None

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
    elif job_type == JobType.REACTIONMINER:
        log.debug(f'Running ReactionMiner job: {job_id}')
        environment = app_config['kubernetes_jobs']['reactionminer']['env']

    elif job_type == JobType.OED_DLKCAT or job_type == JobType.OED_UNIKP or job_type == JobType.OED_CATPRED:
        # Example: "job_info": "{\"input_pairs\":[{\"name\":\"example\",\"sequence\":\"MEDIPDTSRPPLKYVK...\",\"type\":\"FASTA\",\"smiles\":\"OC1=CC=C(C[C@@H](C(O)=O)N)C=C1\"}]}"
        log.debug(f'Running OpenEnzymeDB job: {job_type} - {job_id}')
        log.debug(f'    job_info: {job_info}')
        job_config = json.loads(job_info.replace('\\"', '"'))
        log.debug(f'    job_config: {job_config}')
        job_config_str = json.dumps(job_config['input_pairs']).replace('"', '\\"')
        environment = [{'name': 'OED_INPUT_PAIRS', 'value': job_info.replace('"', '\\"')}]
        log.debug(f'    environment: {environment}')

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

    elif job_type == JobType.NOVOSTOIC_ENZRANK:
        if service.ensure_bucket_exists(job_type):
            upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
            if not upload_result:
                raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
        command = app_config['kubernetes_jobs'][job_type]['command']

    elif job_type == JobType.NOVOSTOIC_DGPREDICTOR:
        if service.ensure_bucket_exists(job_type):
            upload_result = service.upload_file(job_type, f"/{job_id}/in/input.json", job_info.replace('\"', '"').encode('utf-8'))
            if not upload_result:
                raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
        command = app_config['kubernetes_jobs'][job_type]['command']

    elif job_type == JobType.CLEAN:
        # Build up input.FASTA from user input
        job_config = json.loads(job_info.replace('\"', '"'))
        command = CleanService.build_clean_job_command(job_id=job_id, job_info=job_config)
    elif job_type == JobType.MOLLI:
        # Pass path to CORES/SUBS files into the container
        command = app_config['kubernetes_jobs'][job_type]['command']
        job_config = json.loads(job_info.replace('\"', '"'))
        environment = MolliService.build_molli_job_environment(job_id=job_id, job_info=job_config)

    elif job_type == JobType.OED_CHEMINFO:
        # Pass path to CORES/SUBS files into the container
        if service.ensure_bucket_exists(job_type):
            upload_result = service.upload_file(job_type, f"/{job_id}/in/job.json", job_info.replace('\"', '"').encode('utf-8'))
            if not upload_result:
                raise HTTPException(status_code=400, detail="Failed to upload file to MinIO")
        command = app_config['kubernetes_jobs'][job_type]['command']

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

        # Record the parent on the row itself. The relationship was already being
        # passed to the container as an envvar, but nothing persisted it, so a
        # subjob's DB row gave no indication it belonged to anything.
        parent_job_id = job_config['parent_job_id']

        # Pass parent / subjob_ids along as envvars
        environment += app_config['kubernetes_jobs'][job_type]['env']
        environment += [
            { "name": "PARENT_JOB_ID",  "value": job_config['parent_job_id'] },
            { "name": "EZSPEC_UNIDOCK_JOB_ID", "value": job_config['ezspec_unidock_job_id'] },
            { "name": "EZSPEC_INFERENCE_JOB_ID", "value": job_config['ezspec_inference_job_id'] }
        ]

    return PreparedJob(
        job_id=job_id,
        job_info=job_info,
        image_name=image_name,
        command=command,
        environment=environment,
        parent_job_id=parent_job_id,
    )
