# Start from Python 3.10 in the /code directory
FROM python:3.10
WORKDIR /code

# Install PIP dependencies
COPY ./requirements.txt /code/requirements.txt
RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

# Copy in the source code
COPY ./app /code/app

# Copy in versions of the database schema
COPY app/migrations /code/migrations
COPY ./alembic.ini /code/alembic.ini

# Override default command in the container
WORKDIR /code/app
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8080"]

