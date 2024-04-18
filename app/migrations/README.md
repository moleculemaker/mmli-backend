# Database Migrations
Migrations are handled using [Alembic](https://alembic.sqlalchemy.org/en/latest/)

Generic single-database configuration with an async dbapi.

## Create a new Revision:
For example, to create a new revision that will add `new_field` as a column to the `job` table:
```bash
% alembic revision --autogenerate -m "add new_field"
INFO  [alembic.runtime.migration] Context impl PostgresqlImpl.
INFO  [alembic.runtime.migration] Will assume transactional DDL.
INFO  [alembic.ddl.postgresql] Detected sequence named 'job_id_seq' as owned by integer column 'job(id)', assuming SERIAL and omitting
INFO  [alembic.autogenerate.compare] Detected added column 'job.new_field'
  Generating /Users/lambert8/workspace/mmli/mmli-backend/migrations/versions/6a32f8fc32ca_add_new_field.py ...  done
```

You can include `-m comment` with a comment describing the schema change

## Upgrade to Revision
Upgrade to the latest version using `head`, or pass a specific hash instead:
```bash
% alembic upgrade head
INFO  [alembic.runtime.migration] Context impl PostgresqlImpl.
INFO  [alembic.runtime.migration] Will assume transactional DDL.
INFO  [alembic.runtime.migration] Running upgrade fd23b24f6f9a -> 6a32f8fc32ca, add new_field
```

In the database, we see our `new_field` column has been added:
```sql
postgres=# \d job;
                                    Table "public.job"
    Column    |       Type        | Collation | Nullable |             Default             
--------------+-------------------+-----------+----------+---------------------------------
 job_info     | character varying |           |          | 
 email        | character varying |           |          | 
 id           | integer           |           | not null | nextval('job_id_seq'::regclass)
 job_id       | character varying |           |          | 
 run_id       | character varying |           |          | 
 phase        | character varying |           | not null | 
 type         | character varying |           | not null | 
 image        | character varying |           | not null | 
 command      | character varying |           |          | 
 time_created | integer           |           | not null | 
 time_start   | integer           |           | not null | 
 time_end     | integer           |           | not null | 
 deleted      | integer           |           | not null | 
 user_agent   | character varying |           | not null | 
 new_field    | character varying |           |          | 
 ```

## Downgrade to Revision
Downgrade to a specific revision using its hash:
```bash
% alembic downgrade fd23b24f6f9a
INFO  [alembic.runtime.migration] Context impl PostgresqlImpl.
INFO  [alembic.runtime.migration] Will assume transactional DDL.
INFO  [alembic.runtime.migration] Running downgrade 6a32f8fc32ca -> fd23b24f6f9a, add new_field
```

In the database, we see our `new_field` column has been removed:
```sql
postgres=# \d job;
                                    Table "public.job"
    Column    |       Type        | Collation | Nullable |             Default             
--------------+-------------------+-----------+----------+---------------------------------
 job_info     | character varying |           |          | 
 email        | character varying |           |          | 
 id           | integer           |           | not null | nextval('job_id_seq'::regclass)
 job_id       | character varying |           |          | 
 run_id       | character varying |           |          | 
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
