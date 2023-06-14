import io
import os
from minio import Minio
from minio.error import S3Error
from dotenv import load_dotenv

class MinIOService:
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
