# mmli-backend
Unified FastAPI based backend for ChemScraper, (CLEAN job-manager and Molli - future scope)

## The `/v1` API

A self-describing API for running these tools from a script, a notebook, or an agent.
Interactive documentation is at `/v1/docs`; the OpenAPI 3.1 document is at
`/v1/openapi.json`.

There is no authentication. A job is reachable by anyone holding its `job_id`, which is
a server-assigned UUIDv4 and is never listed anywhere. **Treat a `job_id` as a secret.**
`/v1/service-info` states this so a client does not have to guess.

### Find a tool and read what it wants

```bash
BASE=https://mmli.fastapi.mmli1.ncsa.illinois.edu/v1

curl -s $BASE/tools | jq '.itemListElement[] | {identifier, abstract}'
curl -s $BASE/tools/novostoic-optstoic/input-schema | jq .
```

### Run it

```bash
JOB=$(curl -s -X POST $BASE/tools/novostoic-optstoic/jobs \
  -H 'Content-Type: application/json' \
  -H "Idempotency-Key: $(uuidgen)" \
  -d '{"primary_precursor": "MNXM1137670", "target_molecule": "MNXM26"}')

echo "$JOB" | jq -r .job_id
```

`Idempotency-Key` is optional but worth sending: if the request times out, replaying the
same key returns the original job rather than starting a second one.

Inputs are validated against the published schema before anything runs. A rejection
names the offending field:

```json
{
  "type": ".../v1/problems/invalid-input",
  "title": "Input failed schema validation",
  "status": 422,
  "errors": [
    {"pointer": "/primary_precursor", "detail": "5 is not of type 'string'"}
  ]
}
```

### Wait for it, then read the results

```bash
ID=$(echo "$JOB" | jq -r .job_id)

# 409 while running, 200 when finished. Status codes carry the meaning, so there is
# nothing to parse in the polling loop.
until curl -sf -o results.json "$BASE/jobs/$ID/results"; do sleep 10; done
jq . results.json
```

`GET /v1/jobs/$ID` returns status, ISO-8601 timestamps, provenance (including the image
digest that actually ran, once the pod reports it), and a `links` object — so a client
never has to build a URL.

Raw output files are listed at `/v1/jobs/$ID/artifacts`; add `?include=logs` for the
tool's stdout/stderr. `POST /v1/jobs/$ID/cancel` stops a running job.

### Tools that need files

Send one multipart request. The server assigns the `job_id` and stores the files itself:

```bash
curl -X POST $BASE/tools/molli/jobs \
  -F 'inputs={"CORES_FILE_NAME":"cores.cdxml","SUBS_FILE_NAME":"subs.cdxml"};type=application/json' \
  -F 'files=@cores.cdxml' \
  -F 'files=@subs.cdxml'
```

### Errors

Every `/v1` error is [RFC 9457](https://www.rfc-editor.org/rfc/rfc9457) `problem+json`
with a stable `type` URI, so a client can branch on the kind of failure without matching
on prose.

### How this differs from the legacy API

The endpoints under `/{job_type}/...` are unchanged and remain supported. `/v1` differs
in that job ids are server-assigned, inputs are schema-validated, results return 409
rather than `200` with a null body while a job is running, responses carry links, and
errors are problem documents.

## Running the tests

The suite runs inside the same base image the service ships from, so it exercises the
pinned interpreter and library versions rather than whatever is on your machine:

```bash
docker build -f Dockerfile.test -t mmli-backend-test .
docker run --rm mmli-backend-test pytest -q
```

To iterate on tests without rebuilding, mount them in:

```bash
docker run --rm -v "$PWD/tests:/code/tests" mmli-backend-test pytest -q
```

The image is pinned to `linux/amd64`. This is not optional: `python-terrier` depends on
`pytrec-eval-terrier`, which publishes no `aarch64` wheel and cannot build from source,
so an arm64 build fails to install the application's dependencies at all.

### Running the tests against Postgres

By default the suite uses a throwaway SQLite file, so it needs no running services. Set
`TEST_DATABASE_URL` to run the identical tests against the driver production uses:

```bash
docker network create mmli-test-net
docker run -d --name mmli-pg --network mmli-test-net \
  -e POSTGRES_PASSWORD=pw -e POSTGRES_USER=postgres -e POSTGRES_DB=mmli postgres:15

docker run --rm --network mmli-test-net \
  -e TEST_DATABASE_URL="postgresql+asyncpg://postgres:pw@mmli-pg:5432/mmli" \
  mmli-backend-test pytest -q
```

Worth doing whenever the ORM layer changes. SQLite exercises neither `asyncpg` nor
Postgres type handling, so a whole class of driver-level problem is invisible to the
default run.

### Known limitation: the test schema is not the deployed schema

The suite builds its tables with `SQLModel.metadata.create_all()`. Deployments build
them with `alembic upgrade head`. **These do not produce the same schema.** The clearest
example is `job.type`: the migrations declare it as `AutoString` (a plain `varchar`),
while the SQLModel metadata maps the `JobType` enum to a *native Postgres enum*. Under
the migration-built schema an unrecognized value simply matches no rows; under the
metadata-built one the driver raises `InvalidTextRepresentationError`.

So a test can pass or fail for reasons that do not apply in production, in either
direction. Treat behavior that depends on column-level type enforcement as unverified
until it has been checked against a migrated database.

Fixing this properly means building the test schema from migrations, which is not
currently possible on SQLite — several migrations use `ALTER COLUMN ... SET NOT NULL`,
which SQLite does not support — so it would make Postgres mandatory for running tests.
That trade-off has not been made.

### About `tests/characterization/`

These tests pin the **current** behavior of the legacy API, including behavior that is
wrong. Assertions covering known defects carry a `DEFECT:` comment naming the problem.
Their purpose is to make behavioral change visible: when a later change alters one of
these responses, the test diff is the record of that decision. A failure here means
"something changed" — decide whether the change was intended before editing the test.

## ⭐️ Recommended local development (Docker)

### (1/4) Create a `.env`
Create a `.env` from the `env.tpl` file in this repo. The default env is fine without modifications for testing. Change the passwords for production use.
```bash
cp .env.tpl .env
```

### (2/4) Setup a K8 cluster, here we use Minikube 
1. [Install Minikube](https://minikube.sigs.k8s.io/docs/start/?arch=%2Fmacos%2Farm64%2Fstable%2Fbinary+download).
2. Start minikube w/ an external network (defined in this repo's `docker-compose.yml`)
```bash
minikube start --network=mmli-net --driver=docker --memory=24384
```
3. Ensure it's running: `minikube kubectl cluster-info`
```
Kubernetes control plane is running at https://192.168.49.2:8443
CoreDNS is running at https://192.168.49.2:8443/api/v1/namespaces/kube-system/services/kube-dns:dns/proxy
```

4. Apply the necessary confgurations: 
```bash
# Create `mmli` namespace
minikube kubectl -- create ns mmli

# Apply secret and config
minikube kubectl -- apply -f /home/kastan/ncsa/mmli/mmli-backend/app/cfg/local.secret.yaml -n mmli
minikube kubectl -- apply -f /home/kastan/ncsa/mmli/mmli-backend/app/cfg/local.config.yaml -n mmli

# Create PVC needed by molli jobs
minikube kubectl -- apply -f /home/kastan/ncsa/mmli/mmli-backend/chart/weights.pvc.yaml
```

### (3/4) Run Docker Compose build

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

**Test the service works:** Navigate to [`localhost:8080/docs`](localhost:8080/docs) and you should see the FastAPI Swagger docs.

### (4/4) Initialize the databse

Initialize the Postgres database, this creates the SQL tables.
```bash
docker compose exec -w /code rest alembic upgrade head

# You sould see the logs: 
INFO  [alembic.runtime.migration] Context impl PostgresqlImpl.
INFO  [alembic.runtime.migration] Will assume transactional DDL.
INFO  [alembic.runtime.migration] Running upgrade  -> d775ee615d7b, init
INFO  [alembic.runtime.migration] Running upgrade d775ee615d7b -> 88355d0f323b, added moleculecacheentry for caching molecules, modified job schema, added flaggedmolecule for saving flagged molecules
INFO  [alembic.runtime.migration] Running upgrade 88355d0f323b -> e8569ab45dd1, removed moleculecacheentry
INFO  [alembic.runtime.migration] Running upgrade e8569ab45dd1 -> 30b240622d34, add chemical identifier model and table
```
Finally, verify the tables are created: 

1. exec into the database container, running the `pgsql` command.
```bash
docker exec -it mmli-backend-postgresql psql -U postgres mmli
```
2. Run `\d` command to list tables. 
```bash
psql (15.8 (Debian 15.8-1.pgdg120+1))
Type "help" for help.

mmli=# \d 
                     List of relations
 Schema |            Name            |   Type   |  Owner
--------+----------------------------+----------+----------
 public | alembic_version            | table    | postgres
 public | chemical_identifier        | table    | postgres
 public | chemical_identifier_id_seq | sequence | postgres
 public | flaggedmolecule            | table    | postgres
 public | job                        | table    | postgres
(5 rows)
```
3. Check the jobs table `\d job`:
```bash
mmli=# \d job
                        Table "public.job"
    Column    |       Type        | Collation | Nullable | Default
--------------+-------------------+-----------+----------+---------
 job_info     | character varying |           |          |
 email        | character varying |           |          |
 job_id       | character varying |           | not null |
 run_id       | character varying |           |          |
 phase        | character varying |           | not null |
 type         | character varying |           | not null |
 image        | character varying |           |          |
 command      | character varying |           |          |
 time_created | integer           |           | not null |
 time_start   | integer           |           | not null |
 time_end     | integer           |           | not null |
 deleted      | integer           |           | not null |
 user_agent   | character varying |           | not null |
Indexes:
    "job_id_pk" PRIMARY KEY, btree (job_id)
Referenced by:
    TABLE "flaggedmolecule" CONSTRAINT "flaggedmolecule_job_id_fkey" FOREIGN KEY (job_id) REFERENCES job(job_id)
```

🎉 All done! 🎉 Check the Swagger docs for important commands on [`localhost:8080/docs`](localhost:8080/docs).

### How to monitor running jobs

First, submit a job using Curl, Swagger or Postman. E.g.: 

```bash
curl -X POST https://mmli.kastan.ai/aceretro/jobs \
-H "Content-Type: application/json" \
-d '{
  "job_id": "123",
  "run_id": "123",
  "email": "user@gmail.com",
  "job_info": "{\"nuc\": \"hi\", \"CORES_FILE_NAME\": \"hi\", \"SUBS_FILE_NAME\": \"hi\"}"
}'
```

Monitoring the job:
```bash
# after submitted a job, it should create a pod
minikube kubectl -- get pods -A

# Then get details of the pod, including failures
minikube kubectl -- describe pod mmli-job-molli-123456-j4pwd -n mmli

# get the logs from a pod
minikube kubectl -- logs mmli-job-aceretro-222222222222-n5cl4 -n mmli -c job
```

## Local development Setup (without Docker, not recommended)

### (1/3) Configure Environment
Create a `.env` from the `env.tpl` file in this repo. The default env is fine without modifications for testing. Change the passwords for production use.
```bash
cp .env.tpl .env
```

Setting `DEBUG=true` will enable automatically reload the app when the Python source code changes

### (2/3) Install dependencies
Or, you can use Python + pip if you have them installed locally

To install Dependencies:
```bash
# create a new virtual environment, e.g. for conda `conda create -n mmli-backend python=3.10 -y`
# conda activate mmli-backend
pip install -r requirements.txt
```

This will only run the Python app.

⚠️ You must run `MinIO` and `PostgreSQL` yourself. Set their credentials in the `.env` file.

### (3/3) Initialize the databse

Initialize the Postgres database, this initializes the SQL tables with the "init" migration.
```bash
alembic upgrade head
# You should see these logs:
INFO  [alembic.runtime.migration] Context impl PostgresqlImpl.
INFO  [alembic.runtime.migration] Will assume transactional DDL.
INFO  [alembic.runtime.migration] Running upgrade  -> d775ee615d7b, init
...
```

Finally, verify the tables are created: 

1. Run the `pgsql` command.
```bash
psql mmli
```
2. Run `\d` command to list tables. 
```bash
psql (15.8 (Debian 15.8-1.pgdg120+1))
Type "help" for help.

mmli=# \d 
                     List of relations
 Schema |            Name            |   Type   |  Owner
--------+----------------------------+----------+----------
 public | alembic_version            | table    | postgres
 public | chemical_identifier        | table    | postgres
 public | chemical_identifier_id_seq | sequence | postgres
 public | flaggedmolecule            | table    | postgres
 public | job                        | table    | postgres
(5 rows)
```
3. Check the `job` table `\d job`:
```bash
mmli=# \d job
                        Table "public.job"
    Column    |       Type        | Collation | Nullable | Default
--------------+-------------------+-----------+----------+---------
 job_info     | character varying |           |          |
 email        | character varying |           |          |
 job_id       | character varying |           | not null |
 run_id       | character varying |           |          |
 phase        | character varying |           | not null |
 type         | character varying |           | not null |
 image        | character varying |           |          |
 command      | character varying |           |          |
 time_created | integer           |           | not null |
 time_start   | integer           |           | not null |
 time_end     | integer           |           | not null |
 deleted      | integer           |           | not null |
 user_agent   | character varying |           | not null |
Indexes:
    "job_id_pk" PRIMARY KEY, btree (job_id)
Referenced by:
    TABLE "flaggedmolecule" CONSTRAINT "flaggedmolecule_job_id_fkey" FOREIGN KEY (job_id) REFERENCES job(job_id)
```

🎉 All done! 🎉 Check the Swagger docs for important commands on [`localhost:8080/docs`](localhost:8080/docs).

# Database Migrations
Any time that you add, modify, or remove anything in the Job or JobBase classes, this will affect the database schema.

Migrations are handled using [Alembic](https://alembic.sqlalchemy.org/en/latest/)

You can use Alembic to automatically generate a script that will migrate the database to a new schema version.

See the [migrations](migrations/README.md) README for more info
