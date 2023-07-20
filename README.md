# mmli-backend
Unified FastAPI based backend for ChemScraper, (CLEAN job-manager and Molli - future scope)

## Getting Started
You'll need either Docker or Python + pip installed

### Configure Environment
Create a `.env` from the `env.tpl` file in this repo.

You can edit the `.env` file to configure the Python application

For example, setting `DEBUG=true` will enable automatically reload the app when the Python source code changes

### Running in Docker
If you have Docker installed:
```bash
docker compose up -d --build
```

This will run MinIO + PostgreSQL + the Python app

### Running Locally
Or, if you have Python + pip installed locally:
```bash
pip install -r requirements.txt
```

This will only run the Python app.

You must run MinIO + PostgreSQL yourself in this case


## Database Migrations
Migrations are handled using [Alembic](https://alembic.sqlalchemy.org/en/latest/)

See the [migrations](./migrations/README.md) README for more info
