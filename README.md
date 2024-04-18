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
Or, you can use Python + pip if you have them installed locally

To install Dependencies:
```bash
pip install -r requirements.txt
```


This will only run the Python app.

You must run MinIO + PostgreSQL yourself in this case


## Database Migrations
Any time that you add, modify, or remove anything in the Job or JobBase classes, this will affect the database schema.

Migrations are handled using [Alembic](https://alembic.sqlalchemy.org/en/latest/)

You can use Alembic to automatically generate a script that will migrate the database to a new schema version.

See the [migrations](app/migrations/README.md) README for more info


### Initialization
When running the first time, you'll need to create the needed database tables.

In Docker:
```bash
docker compose exec -w /code -it rest alembic upgrade head
```

Outside of Docker:
```bash
alembic upgrade head
```

This will apply the "init" migration and create the needed database tables:
```bash
INFO  [alembic.runtime.migration] Context impl PostgresqlImpl.
INFO  [alembic.runtime.migration] Will assume transactional DDL.
INFO  [alembic.runtime.migration] Running upgrade  -> d775ee615d7b, init
```

In PostgreSQL, you should see the new table has been created:
```bash
% docker exec -it mmli-backend-postgresql psql -U postgres mmli
psql (15.3 (Debian 15.3-1.pgdg120+1))
Type "help" for help.

mmli=# \d
               List of relations
 Schema |      Name       |   Type   |  Owner   
--------+-----------------+----------+----------
 public | alembic_version | table    | postgres
 public | job             | table    | postgres
 public | job_id_seq      | sequence | postgres
(3 rows)

mmli=# \d job;
                                    Table "public.job"
    Column    |       Type        | Collation | Nullable |             Default             
--------------+-------------------+-----------+----------+---------------------------------
 job_info     | character varying |           |          | 
 email        | character varying |           |          | 
 job_id       | character varying |           |          | 
 run_id       | character varying |           |          | 
 id           | integer           |           | not null | nextval('job_id_seq'::regclass)
 phase        | character varying |           | not null | 
 type         | character varying |           | not null | 
 image        | character varying |           | not null | 
 command      | character varying |           |          | 
 time_created | integer           |           | not null | 
 time_start   | integer           |           | not null | 
 time_end     | integer           |           | not null | 
 deleted      | integer           |           | not null | 
 user_agent   | character varying |           | not null | 
Indexes:
    "job_pkey" PRIMARY KEY, btree (id)
```
