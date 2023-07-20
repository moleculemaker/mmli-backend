# Start from Python 3.10 in the /code directory
FROM python:3.10
WORKDIR /code

# Install PIP dependencies
COPY ./requirements.txt /code/requirements.txt
RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

# Copy in the source code
COPY ./app /code/app
COPY ./migrations /code/migrations

# Set URLS to postgresql / minio (these can be overridden at runtime)
ENV DATABASE_URL=postgresql://postgres:postgres@localhost:5432/mmli

# Override default command in the container
WORKDIR /code/app
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8080"]

