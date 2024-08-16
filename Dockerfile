# Start from Python 3.10 in the /code directory
FROM condaforge/mambaforge:24.1.2-0
WORKDIR /code

# Install PIP dependencies
COPY ./requirements.txt /code/requirements.txt
RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt
RUN mamba install openbabel

# Copy in the source code
COPY ./app /code/app

# Copy in versions of the database schema
COPY migrations /code/migrations
COPY ./alembic.ini /code/alembic.ini

# Override default command in the container
WORKDIR /code/app
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8080"]