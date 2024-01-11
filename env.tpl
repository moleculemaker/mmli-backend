# App will automatically reload on change when DEBUG=True
DEBUG="true"

# Set ChemScraper backend API URL
CHEMSCRAPER_API_BASE_URL=https://chemscraper.backend.staging.mmli1.ncsa.illinois.edu

# Set ChemScraper backend API URL
CHEMSCRAPER_FRONTEND_URL=http://localhost:4200

# Email server config
EMAIL_SERVER=smtp.ncsa.uiuc.edu
EMAIL_FROM_EMAIL=devnull+chemscraper@ncsa.illinois.edu
EMAIL_FROM_NAME=no-reply-ChemScraper


# Set some default values (override below if needed)
POSTGRES_USER=postgres
POSTGRES_PASSWORD=postgres
POSTGRES_DB=mmli

# Set MinIO default values (override below if needed)
MINIO_ROOT_USER=minioadmin
MINIO_ROOT_PASSWORD=minioadmin
MINIO_ACCESS_KEY=${MINIO_ROOT_USER}
MINIO_SECRET_KEY=${MINIO_ROOT_PASSWORD}

# For local development with your own running instances (without Docker)
#MINIO_SERVER=localhost:9000
#POSTGRES_SERVER=localhost:5432

# If your db/minio are running in Docker, then localhost within the container does not resolve to the right service
# Use "host.docker.internal" as hostname in this case instead
MINIO_SERVER=host.docker.internal:9000
POSTGRES_SERVER=host.docker.internal:5432

# Alternatively, if you are using Docker compose for local deployment,
# then hostname in the connection URL can be left as "postgresql" or "minio"
#POSTGRES_SERVER=postgresql:5432
#MINIO_SERVER=minio:9000

# If you are using Kubernetes deployment, and you have already deployed the Postgresql service in a Pod, then
# you can use the DNS as postgresql+asyncpg://<POSTGRES_USER>:<POSTGRES_PASSWORD>@SERVICE_NAME.NAMESPACE.svc.cluster.local:5432/<POSTGRES_DATABASE>
#POSTGRES_SERVER=SERVICE_NAME.NAMESPACE.svc.cluster.local:5432
#MINIO_SERVER=SERVICE_NAME.NAMESPACE.svc.cluster.local:9000

# Set database URL based on other envs
SQLALCHEMY_DATABASE_URL=postgresql+asyncpg://${POSTGRES_USER}:${POSTGRES_PASSWORD}@${POSTGRES_SERVER}/${POSTGRES_DB}
