import io
import logging

from minio import Minio
from minio.error import S3Error
from urllib.parse import urlparse, urlunparse

from config import app_config, get_logger, MINIO_SERVER, MINIO_ACCESS_KEY, MINIO_SECRET_KEY

log = get_logger(__name__)

class MinIOService:
    minio_api_baseURL = app_config['minio']['apiBaseUrl']

    def __init__(self):
        log.info(f'Setting up MinIO client: {MINIO_SERVER} {MINIO_ACCESS_KEY}:*********')
        self.client = Minio(
            MINIO_SERVER,
            access_key=MINIO_ACCESS_KEY,
            secret_key=MINIO_SECRET_KEY,
            secure=False
        )

    def get_file(self, bucket_name, object_name):
        try:
            log.info(f"Fetching file:   bucket_name={bucket_name}   object_name={object_name}")
            data = self.client.get_object(bucket_name, object_name)
            return data.read()
        except S3Error as err:
            log.error("Error: ", err)

    def upload_file(self, bucket_name, object_name, file_content):
        if isinstance(file_content, str):
            file_content_bytes = file_content.encode('utf-8')
            file_content = io.BytesIO(file_content_bytes)

        try:
            # Upload file to the specified bucket
            result = self.client.put_object(
                bucket_name, object_name, io.BytesIO(file_content), len(file_content)
            )
            log.info(
                "File uploaded successfully. Object name: {}, Etag: {}".format(
                    object_name, result.etag
                )
            )
            return True
        except S3Error as err:
            log.error("Error: ", err)
            return False

    def upload_file_from_str(self, bucket_name, object_name, contents):
        value_as_bytes = contents.encode('utf-8')
        value_as_a_stream = io.BytesIO(value_as_bytes)
        self.upload_file()

    def list_files(self, bucket_name, path, recursive=False):
        return self.client.list_objects(bucket_name, prefix=path, recursive=recursive)
        
    def get_file_urls(self, bucket_name, path):
        try:
            urls = []
            objects = self.client.list_objects(bucket_name, prefix=path, recursive=True)
            for obj in objects:
                url = self.client.presigned_get_object(bucket_name, obj.object_name)
                url = url.split('?', 1)[0]
                parsed_url = urlparse(url)
                minio_api_url = urlunparse(
                    ("https", self.minio_api_baseURL, parsed_url.path, parsed_url.params, parsed_url.query, parsed_url.fragment)
                )
                # urls.append(url)
                urls.append(minio_api_url)
            return urls
        except S3Error as err:
            log.error("Error: ", err)

    def ensure_bucket_exists(self, bucket_name):
        if self.client.bucket_exists(bucket_name):
            log.debug(f"{bucket_name} already exists. Using existing bucket: {bucket_name}")
            return True
        else:
            log.debug(f"{bucket_name} does not exist. Creating new bucket: {bucket_name}")
            try:
                self.client.make_bucket(bucket_name=bucket_name)
                return True
            except:
                return False

    def check_file_exists(self, bucket_name, object_name):
        try:
            self.client.stat_object(bucket_name, object_name)
            return True
        except S3Error as err:
            if err.code == "NoSuchKey":
                return False
            else:
                log.error("Error: ", err)
                return False
