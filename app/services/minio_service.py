import io
import os
from minio import Minio
from minio.error import S3Error
from dotenv import load_dotenv
from urllib.parse import urlparse, urlunparse

class MinIOService:
    minio_api_baseURL = "minioapi.mmli.fastapi.staging.mmli1.ncsa.illinois.edu"

    def __init__(self):
        load_dotenv()
        self.client = Minio(
            os.environ.get("MINIO_SERVER"),
            access_key=os.environ.get("MINIO_ACCESS_KEY"),
            secret_key=os.environ.get("MINIO_SECRET_KEY"),
            secure=False
        )

    def get_file(self, bucket_name, object_name):
        try:
            data = self.client.get_object(bucket_name, object_name)
            return data.read()
        except S3Error as err:
            print("Error: ", err)

    def upload_file(self, bucket_name, object_name, file_content):
        try:
            # Upload file to the specified bucket
            result = self.client.put_object(
                bucket_name, object_name, io.BytesIO(file_content), len(file_content)
            )
            print(
                "File uploaded successfully. Object name: {}, Etag: {}".format(
                    object_name, result.etag
                )
            )
            return True
        except S3Error as err:
            print("Error: ", err)
            return False
        
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
            print("Error: ", err)

    def check_file_exists(self, bucket_name, object_name):
        try:
            self.client.stat_object(bucket_name, object_name)
            return True
        except S3Error as err:
            if err.code == "NoSuchKey":
                return False
            else:
                print("Error: ", err)
                return False
