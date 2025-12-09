from minio import Minio


if __name__ == "__main__":
  client = Minio(
    'minio.mmli.fastapi.mmli1.ncsa.illinois.edu:9000',
    access_key='8aa3c114d36ba20212bc307e2f321c6b',
    secret_key='53a450ba3a2c72b876bccfcfafabf16d',
    secure=False
  )
  
  bucket_name = "oed-dlkcat"
  path = "/b391d8ee04ef4180aea3367e66e27757/out/"
  recursive = False
  
  files = client.list_objects(bucket_name, prefix=path, recursive=recursive)
  print(list(files))
  # print(files)
  # print(files)
    