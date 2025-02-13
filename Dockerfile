# Start from Python 3.10 in the /code directory
FROM --platform=linux/amd64 condaforge/mambaforge:24.1.2-0
WORKDIR /code

# Set noninteractive installation and timezone
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Install build dependencies and OpenJDK
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    zlib1g-dev \
    openjdk-11-jdk \
    tzdata \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME environment variable
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-arm64
ENV PATH=$PATH:$JAVA_HOME/bin

# Install conda dependencies first
RUN mamba install -y openbabel

# Install PIP dependencies
COPY ./requirements.txt /code/requirements.txt
RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

# Copy in the source code
COPY ./app /code/app

# Copy in versions of the database schema
COPY migrations /code/migrations
COPY ./alembic.ini /code/alembic.ini

# Override default command in the container
WORKDIR /code/app
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8080"]