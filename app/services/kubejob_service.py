import logging
import sys

#import global_vars
#from global_vars import config, log
from kubernetes import watch, client, config as kubeconfig
from kubernetes.client.rest import ApiException
import os
import yaml
import json
from jinja2 import Template
import uuid

import time
import threading

from requests import HTTPError
from sqlmodel import Session
from sqlmodel.ext.asyncio.session import AsyncSession

import config
from config import app_config
from models.enums import JobStatus
from models.sqlmodel import db
from models.sqlmodel.db import get_session, engine
from models.sqlmodel.models import Job

#from dbconnector import DbConnector

log = config.get_logger(__name__)

## Load Kubernetes cluster config. Unhandled exception if not in Kubernetes environment.
try:
    kubeconfig.load_incluster_config()
    log.info('Successfully loaded in-cluster config!')
except Exception as e1:
    try:
        config_file_path = app_config['kubernetes_jobs']['defaults']['kubeconfig']
        log.info('In-cluster config failed: Falling back to provided kubeconfig path: ' + config_file_path)
        kubeconfig.load_kube_config(config_file=config_file_path)
    except Exception as e2:
        log.fatal('Failed to get any cluster config: ', e2)
        sys.exit(1)
#configuration = client.Configuration()
#api_batch_v1 = client.BatchV1Api(client.ApiClient(configuration))
#api_v1 = client.CoreV1Api(client.ApiClient(configuration))
api_batch_v1 = client.BatchV1Api()
api_v1 = client.CoreV1Api()



# watched_namespaces = ["test"]
# for namespace in watched_namespaces:
#    count = 10
#    w = watch.Watch()
#    for event in w.stream(custom.list_namespaced_custom_object,
#                          group=CRD_GROUP_NAME,
#                          version=CRD_VERSION_V1,
#                          namespace=namespace,
#                          plural=CRD_USERAPPS_PLURAL,
#                          _request_timeout=60):
#        print("Event: %s %s" % (event['type'], event['object']['metadata']['name']))
#        count -= 1
#        time.sleep(10)
#        if not count:
#            w.stop()

async def get_db() -> AsyncSession:
    async with get_session() as db:
        try:
            yield db
        except:
            # if we fail somehow rollback the connection
            log.warn("We somehow failed in a DB operation and auto-rollbacking...")
            await db.rollback()
            raise

class KubeEventWatcher:

    def __init__(self):
        self.logger = log
        #self.logger.setLevel('DEBUG')
        self.thread = threading.Thread(target=self.run, name='kube-event-watcher', daemon=True)
        # Get global instance of the job handler database interface
        self.db = db.create_db_engine()

        self.stream = None
        self.logger.info('Starting KubeWatcher')
        self.thread.start()
        self.logger.info('Started KubeWatcher')

    def run(self):
        # Ignore kube-system namespace
        # TODO: Parameterize this?
        ignored_namespaces = ['kube-system']
        self.logger.info('KubeWatcher watching all namespaces except for: ' + str(ignored_namespaces))

        # TODO: Parameterize this?
        required_labels = {
            'type': 'uws-job'
        }
        self.logger.info('KubeWatcher looking for required labels: ' + str(required_labels))

        timeout_seconds = 0
        resource_version = ''
        k8s_event_stream = None

        w = watch.Watch()

        while True:
            time.sleep(1)
            self.logger.info('KubeWatcher is connecting...')
            try:
                # List all pods in watched namespace to get resource_version
                namespaced_jobs: client.V1JobList = api_batch_v1.list_namespaced_job(namespace=get_namespace())
                resource_version = namespaced_jobs.metadata.resource_version if namespaced_jobs.metadata.resource_version else resource_version

                # Then, watch for new events using the most recent resource_version
                # Resource version is used to keep track of stream progress (in case of resume/retry)
                k8s_event_stream = w.stream(func=api_batch_v1.list_namespaced_job,
                                            namespace=get_namespace(),
                                            resource_version=resource_version,
                                            timeout_seconds=timeout_seconds)

                self.logger.info('KubeWatcher connected!')

                # Parse events in the stream for Pod phase updates
                for event in k8s_event_stream:
                    resource_version = event['object'].metadata.resource_version

                    # Skip Pods in ignored namespaces
                    if event['object'].metadata.namespace in ignored_namespaces:
                        self.logger.debug('Skipping event in excluded namespace')
                        continue

                    # Examine labels, ignore if not uws-job
                    # self.logger.debug('Event recv\'d: %s' % event)
                    labels = event['object'].metadata.labels

                    if labels is None and len(required_labels) > 0:
                        self.logger.warning(
                            'WARNING: Skipping due to missing label(s): ' + str(required_labels))
                        continue

                    missing_labels = [x for x in required_labels if x not in labels]
                    if len(missing_labels) > 0:
                        self.logger.warning(
                            'WARNING: Skipping due to missing required label(s): ' + str(missing_labels))
                        continue

                    # TODO: lookup associated userapp using resource name
                    name = event['object'].metadata.name

                    # Parse name into userapp_id + ssid
                    segments = name.split('-')
                    if len(segments) < 3:
                        self.logger.warning('WARNING: Invalid number of segments -  JobName=%s' % name)
                        continue

                    # uws-job-jobid => we want last segment
                    job_id = segments[-1]

                    type = event['type']
                    status = event['object'].status
                    conditions = status.conditions

                    # Calculate new status
                    self.logger.debug(f'Event: job_id={job_id}   type={type}   status={status}')
                    new_phase = None
                    if conditions is None:
                        new_phase = 'processing'
                    elif len(conditions) > 0 and conditions[0].type == 'Complete':
                        new_phase = 'completed'
                    elif status.failed > 0:
                        new_phase = 'error'
                    else:
                        self.logger.info(f'>> Skipped job update: {job_id}-> {new_phase}')
                        self.logger.debug(f'>> Status: {str(status)}')

                    # Write status update back to database
                    if new_phase is not None:
                        self.logger.debug('Updating job phase: %s -> %s' % (job_id, new_phase))

                        with Session(engine) as session:
                            db.update_job_phase(session, job_id, JobStatus(new_phase))
                            self.logger.info('Updated job phase: %s -> %s' % (job_id, new_phase))

            except (ApiException, HTTPError) as e:
                self.logger.error('HTTPError encountered - KubeWatcher reconnecting to Kube API: %s' % str(e))
                if k8s_event_stream:
                    k8s_event_stream.close()
                k8s_event_stream = None
                if e.status == 410:
                    # Resource too old
                    resource_version = ''
                    self.logger.warning("Resource too old (410) - reconnecting: " + str(e))
                time.sleep(2)
                continue
            except Exception as e:
                self.logger.error('Unknown exception - KubeWatcher reconnecting to Kube API: %s' % str(e))
                if k8s_event_stream:
                    k8s_event_stream.close()
                k8s_event_stream = None
                time.sleep(2)
                continue

    def is_alive(self):
        return self.thread and self.thread.is_alive()

    def close(self):
        if self.stream:
            self.stream.close()
            self.stream = None
        if self.is_alive():
            self.thread.join(timeout=3)
            self.thread = None


def get_namespace():
    # When running in a pod, the namespace should be determined automatically,
    # otherwise we assume the local development is in the default namespace
    try:
        with open('/var/run/secrets/kubernetes.io/serviceaccount/namespace', 'r') as file:
            namespace = file.read().replace('\n', '')
    except:
        namespace = app_config['kubernetes_jobs']['defaults']['namespace']

    return namespace


def generate_uuid():
    return str(uuid.uuid4()).replace("-", "")


def get_job_name_from_id(job_type, job_id):
    return 'mmli-job-{}-{}'.format(job_type, job_id)


def get_job_root_dir_from_id(job_id):
    return os.path.join(app_config['kubernetes_jobs']['defaults']['workingVolume']['mountPath'], 'jobs', job_id)


def list_job_output_files(job_id):
    job_filepaths = []
    try:
        job_output_dir = os.path.join(get_job_root_dir_from_id(job_id), 'out')
        log.debug(f'Listing job files ({job_id}) [{job_output_dir}]:')
        if os.path.isdir(job_output_dir):
            for dirpath, dirnames, filenames in os.walk(job_output_dir):
                for filename in filenames:
                    filepath = os.path.join(dirpath, filename)
                    if not filename.startswith('.') and os.path.isfile(filepath):
                        job_filepaths.append(filepath)
                        log.debug(filepath)
    except Exception as e:
        log.error(str(e))
        raise e
    return job_filepaths


def list_jobs(job_id=None):
    jobs = []
    response = {
        'jobs': jobs,
        'status': config.STATUS_OK,
        'message': '',
    }
    try:
        namespace = get_namespace()
        if job_id:
            api_response = api_batch_v1.list_namespaced_job(
                namespace=namespace,
                label_selector=f'jobId={job_id}',
            )
        else:
            api_response = api_batch_v1.list_namespaced_job(
                namespace=namespace,
                label_selector=f'type=mmli-job',
            )
        ## ref: https://github.com/kubernetes-client/python/blob/549482d8842a47c8b51122422e879e7cc497bf88/kubernetes/docs/V1Job.md
        for item in api_response.items:
            envvars = []
            for envvar in item.spec.template.spec.containers[0].env:
                envvars.append({
                    'name': envvar.name,
                    'value': envvar.value,
                })
            command = []
            for command_element in item.spec.template.spec.containers[0].command:
                command.append(command_element.strip())
            job = {
                'name': item.metadata.name,
                ## CreationTimestamp is a timestamp representing the server time when this object was created.
                'creation_time': item.metadata.creation_timestamp,
                'job_id': item.metadata.labels['jobId'],
                'run_id': item.metadata.labels['runId'],
                'owner_id': item.metadata.labels['ownerId'],
                'command': command,
                'environment': envvars,
                'output_files': list_job_output_files(item.metadata.labels['jobId']),
                'status': {
                    ## active: The number of actively running pods.
                    'active': True if item.status.active else False,
                    ## start_time: Represents time when the job was acknowledged by the job controller.
                    'start_time': item.status.start_time,
                    ## completion_time: Represents time when the job was completed.
                    'completion_time': item.status.completion_time,
                    ## succeeded: The number of pods which reached phase Succeeded.
                    'succeeded': True if item.status.succeeded else False,
                    ## failed: The number of pods which reached phase Failed.
                    'failed': True if item.status.failed else False,
                },
            }
            jobs.append(job)
        response['jobs'] = jobs
    except Exception as e:
        msg = str(e)
        log.error(msg)
        response['status'] = config.STATUS_ERROR
        response['message'] = msg
    return response


def delete_job(job_id: str) -> None:
    namespace = get_namespace()
    job_name = get_job_name_from_id(job_id)
    # config_map_name = f'''{job_name}-positions'''
    body = client.V1DeleteOptions(propagation_policy='Background')
    api_response = None
    try:
        api_response = api_batch_v1.delete_namespaced_job(
            namespace=namespace,
            name=job_name,
            body=body,
        )
    except ApiException as e:
        ## No need to raise exception for 404 not found errors
        log.error(f'''Error deleting Kubernetes Job: {e}''')
        if api_response:
            log.error(json.dumps(api_response))
    except Exception as e:
        log.error(f'''Error deleting Kubernetes Job: {e}''')
        raise
    # try:
    #     api_response = api_v1.delete_namespaced_config_map(
    #         namespace=namespace,
    #         name=config_map_name,
    #         body=body,
    #     )
    # except ApiException as e:
    #     ## No need to raise exception for 404 not found errors
    #     log.error(f'''Error deleting Kubernetes ConfigMap: {e}''')
    #     if api_response:
    #         log.error(json.dumps(api_response))
    # except Exception as e:
    #     log.error(f'''Error deleting Kubernetes ConfigMap: {e}''')
    #     raise


def create_job(job_type, image_name=None, command=None, job_id=None, run_id=None, owner_id=None, replicas=1, environment=[]):
    # TODO: run test/echo job when job_type == 'default'?
    if job_type not in app_config['kubernetes_jobs']:
        log.error(f'Cannot find config for job_type={job_type}')
        return {
            'job_id': None,
            'message': None,
            'command': None,
            'image': None,
            'status': config.STATUS_ERROR,
        }

    image_name = app_config['kubernetes_jobs'][job_type]['image']

    try:
        pullSecrets = app_config['kubernetes_jobs'][job_type]['imagePullSecrets']
        log.debug(f"Using image pullSecrets={pullSecrets}")
    except:
        pullSecrets = None
        log.debug(f"Using default: pullSecrets=None")

    # Set the image pull policy if specified; otherwise set semi-intelligently
    try:
        pullPolicy = app_config['kubernetes_jobs'][job_type]['imagePullPolicy']
    except:
        pullPolicy = 'Always' if image_name in [':latest', ':dev'] else 'IfNotPresent'

    response = {
        'job_id': None,
        'message': None,
        'image': image_name,
        'command': command,
        'status': config.STATUS_OK,
    }

    try:
        namespace = get_namespace()
        if not job_id:
            job_id = generate_uuid()
        if not run_id:
            run_id = job_id
        job_name = get_job_name_from_id(job_type, job_id)
        # config_map_name = f'''{job_name}-positions'''
        job_root_dir = get_job_root_dir_from_id(job_id)
        job_output_dir = os.path.join(job_root_dir, 'out')

        # Set environment-specific configuration for Job definition
        templateFile = "job.tpl.yaml"
        project_subpath = ''

        # all_volumes = [{
        #     'name': 'positions-volume',
        #     'configMapName': config_map_name,
        #     'mountPath': '/etc/cutout/positions.csv',
        #     'subPath': 'positions.csv',
        #     'readOnly': True,
        # }]
        all_volumes = []
        for volume in app_config['kubernetes_jobs']['defaults']['volumes']:
            all_volumes.append(volume)

        jobCompleteApiUrl = f'''{app_config['server']['protocol']}://{os.path.join(
            app_config['server']['hostName'],
            app_config['server']['basePath'],
            app_config['server']['apiBasePath'],
            'uws/report'
        )}'''

        with open(os.path.join(os.path.dirname(__file__), 'templates', templateFile)) as f:
            templateText = f.read()
        template = Template(templateText)


        job_body = yaml.safe_load(template.render(
            name=job_name,
            runId=run_id,
            jobId=job_id,
            ownerId=owner_id,
            namespace=namespace,
            backoffLimit=0,
            replicas=replicas,
            ## CAUTION: Discrepancy between the UID of the image user and the UWS API server UID
            ##          will create permissions problems. For example, if the job UID is 1001 and
            ##          the server UID is 1000, then files created by the job will not in general
            ##          allow the server to delete them when cleaning up deleted jobs.
            image={
                'repository': image_name,
                'pull_policy': pullPolicy,
                'pull_secrets': pullSecrets
            },
            command=command,

            environment=environment,
            uws_root_dir=app_config['kubernetes_jobs']['defaults']['workingVolume']['mountPath'],
            job_output_dir=job_output_dir,
            project_subpath=project_subpath,
            securityContext=app_config['kubernetes_jobs']['defaults']['securityContext'],
            workingVolume=app_config['kubernetes_jobs']['defaults']['workingVolume'],
            volumes=all_volumes,
            resources=app_config['kubernetes_jobs']['defaults']['resources'],
            # apiToken=config['jwt']['hs256Secret'],
            apiToken='dummy',
            jobCompleteApiUrl=jobCompleteApiUrl,
            # releaseName=os.environ.get('RELEASE_NAME', 'dummy'),
            ttlSecondsAfterFinished=app_config['kubernetes_jobs']['defaults']['ttlSecondsAfterFinished'],
            activeDeadlineSeconds=app_config['kubernetes_jobs']['defaults']['activeDeadlineSeconds'],
        ))
        log.debug("Job {}:\n{}".format(job_name, yaml.dump(job_body, indent=2)))
        api_response = api_batch_v1.create_namespaced_job(
            namespace=namespace, body=job_body
        )
        response['job_id'] = job_id

        log.debug(f"Job {job_name} created")
    # TODO: Is there additional information to obtain from the ApiException?
    # except ApiException as e:
    #     msg = str(e)
    #     log.error(msg)
    #     response['status'] = global_vars.STATUS_ERROR
    #     response['message'] = msg
    except ApiException as e:
        msg = str(e)
        log.error(f'Failed creating Kubernetes job: {msg}')
        response['message'] = msg
        response['status'] = config.STATUS_ERROR
    return response