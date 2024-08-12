# mmli-backend
Unified FastAPI based backend for ChemScraper, (CLEAN job-manager and Molli - future scope)

## Recommended local development (Docker)

### (1/3) Create a `.env`
Create a `.env` from the `env.tpl` file in this repo. The default env is fine without modifications for testing. Change the passwords for production use.
```bash
cp .env.tpl .env
```

### (2/3) Setup a K8 cluster, here we use Minikube 
1. [Install Minikube](https://minikube.sigs.k8s.io/docs/start/?arch=%2Fmacos%2Farm64%2Fstable%2Fbinary+download).
2. Follow the instructions to set it up, e.g. `minikube start`
3. Ensure it's running: `minikube kubectl cluster-info`
```
Kubernetes control plane is running at https://192.168.49.2:8443
CoreDNS is running at https://192.168.49.2:8443/api/v1/namespaces/kube-system/services/kube-dns:dns/proxy
```

### (3/3) Run Docker Compose build

Edit the `docker-compose.yml` to expose your kube config. In our case, `minikube` requires 3 values: `ca.crt`, `client.crt`, `client.key`. 

Copy-paste this into the `docker-compose.yml`, under the `rest` container:

> ⚠️ Note: I had problems with `${HOME}` and had to provide full absolute paths manually; e.g. replace `${HOME}` with `/home/username`. ⚠️

```
rest:
    container_name: mmli-backend
    
    ...

    volumes:
        - ./app:/code/app
        - ./migrations:/code/migrations
        - ${HOME}/.kube/config:/opt/kubeconfig
        - ${HOME}/.minikube/ca.crt:/home/kastan/.minikube/ca.crt
        - ${HOME}/.minikube/profiles/minikube/client.crt:${HOME}/.minikube/profiles/minikube/client.crt
        - ${HOME}/.minikube/profiles/minikube/client.key:${HOME}/.minikube/profiles/minikube/client.key
        
```

Finally start the compose. Monitor for errors from `mmli-backend` in the logs.

```bash
docker compose up --build # optionally add -d for detached
```
This will run `MinIO` + `PostgreSQL` + the Python app `mmli-backend`.

**Test the service works:** Navigate to [`localhost:8080/docs`](localhost:8080/docs) and you should see the FastAPI Swagger docs 🎉 Done! 🎉

## Local development Setup without Docker (local install)

### Configure Environment
Create a `.env` from the `env.tpl` file in this repo. The default env is fine without modifications for testing. Change the passwords for production use.
```bash
cp .env.tpl .env
```

Setting `DEBUG=true` will enable automatically reload the app when the Python source code changes

### 
Or, you can use Python + pip if you have them installed locally

To install Dependencies:
```bash
# create a new virtual environment, e.g. for conda `conda create -n mmli-backend python=3.10 -y`
# conda activate mmli-backend
pip install -r requirements.txt
```

This will only run the Python app.

You must run `MinIO` and `PostgreSQL` yourself. Set their credentials in the `.env` file.


# Database Migrations
Any time that you add, modify, or remove anything in the Job or JobBase classes, this will affect the database schema.

Migrations are handled using [Alembic](https://alembic.sqlalchemy.org/en/latest/)

You can use Alembic to automatically generate a script that will migrate the database to a new schema version.

See the [migrations](migrations/README.md) README for more info


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
