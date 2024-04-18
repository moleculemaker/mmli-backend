import os
import logging
import yaml

from collections import ChainMap
from dotenv import load_dotenv


def get_logger(name):
    logger = logging.getLogger(name)
    try:
        logger.setLevel(SERVER_LOGLEVEL.upper())
    except:
        logger.setLevel('INFO')
    return logger


# Load env from .env file, if present
load_dotenv()

# Configuration files for all application settings
CONFIG_FILEPATH = os.getenv('CONFIG_FILEPATH', 'cfg/config.yaml')
SECRET_FILEPATH = os.getenv('SECRET_FILEPATH', 'cfg/secret.yaml')

# Load configuration + secrets files
with open(CONFIG_FILEPATH, "r") as server_yaml_file:
    # app_config takes on the format in cfg/config.yaml
    app_config = yaml.load(server_yaml_file, Loader=yaml.FullLoader)

    with open(SECRET_FILEPATH, "r") as secret_yaml_file:
        # app_secrets takes on the format in cfg/secrets.yaml
        secrets_file = yaml.load(secret_yaml_file, Loader=yaml.FullLoader)

        # app_config is a merged dictionary contains both config_file and secrets_file
        app_config['minio']['accessKey'] = secrets_file['minio']['accessKey']
        app_config['minio']['secretKey'] = secrets_file['minio']['secretKey']
        app_config['auth']['hcaptcha_secret'] = secrets_file['auth']['hcaptcha_secret']

        app_config['db'] = {}
        app_config['db']['url'] = secrets_file['db']['url']

# Configure logging
logging.basicConfig(format='%(asctime)s [%(name)-12s] %(levelname)-8s %(message)s')

log = get_logger(__name__)
log.debug('Server configuration: ', app_config)

# Override app_config with some individual environment variables
# TODO: Is this needed??

SERVER_LOGLEVEL = os.getenv('LOGLEVEL', app_config['server']['loglevel'])

DEBUG = SERVER_LOGLEVEL.lower() == 'debug'
SQLALCHEMY_DATABASE_URL = os.getenv("SQLALCHEMY_DATABASE_URL", app_config['db']['url'])
MINIO_SERVER = os.getenv("MINIO_SERVER", app_config['minio']['server'])
MINIO_ACCESS_KEY = os.getenv("MINIO_ACCESS_KEY", app_config['minio']['accessKey'])
MINIO_SECRET_KEY = os.getenv("MINIO_SECRET_KEY", app_config['minio']['secretKey'])

CHEMSCRAPER_API_BASE_URL = os.getenv('CHEMSCRAPER_API_BASE_URL', app_config['external']['chemscraper']['apiBaseUrl'])

# OAUTH_USERINFO_URL = os.getenv('OAUTH_USERINFO_URL', app_config['auth']['userInfoUrl'])
# OAUTH_COOKIE_NAME = os.getenv('OAUTH_COOKIE_NAME', app_config['auth']['cookieName'])

# TODO: currently unused?
# HCAPTCHA_SECRET = os.getenv('HCAPTCHA_SECRET', app_config['auth']['hcaptcha_secret'])


STATUS_OK = 'ok'
STATUS_ERROR = 'error'

# MMLI_SHARED_MOUNTPATH = os.getenv('MMLI_SHARED_MOUNTPATH', '/mnt/mmli/shared/')
