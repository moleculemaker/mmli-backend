"""Importing config must not write the contents of the secrets file to the log.

config.py logged the whole merged app_config at INFO, so every backend start (and every
job's init and postjob container, which use the same image) printed the MinIO secret key,
the hCaptcha secret and the database URL with its password into the pod log.
"""
import os
import pathlib
import subprocess
import sys

APP_DIR = pathlib.Path(__file__).resolve().parent.parent / "app"

SECRETS = {
    "minio_accessKey": "marker-access-key-3f9a",
    "minio_secretKey": "marker-secret-key-7c21",
    "auth_hcaptcha_secret": "marker-hcaptcha-5b8e",
    "database_url": "postgresql+asyncpg://postgres:marker-db-password-9d4c@db:5432/mmli",
}


def test_importing_config_does_not_log_secrets(tmp_path):
    secret_file = tmp_path / "secret.yaml"
    secret_file.write_text("".join(f'{k}: "{v}"\n' for k, v in SECRETS.items()))

    env = {
        **os.environ,
        "CONFIG_FILEPATH": str(APP_DIR / "cfg" / "config.yaml"),
        "SECRET_FILEPATH": str(secret_file),
    }
    # Keep the test's own overrides from masking what config.py reads from the file.
    for name in ("SQLALCHEMY_DATABASE_URL", "MINIO_ACCESS_KEY", "MINIO_SECRET_KEY"):
        env.pop(name, None)

    result = subprocess.run(
        [sys.executable, "-c", "import config; print(config.MINIO_SECRET_KEY == 'marker-secret-key-7c21')"],
        cwd=APP_DIR, env=env, capture_output=True, text=True, check=True,
    )

    # The secrets were loaded; they just must not appear in the output.
    assert result.stdout.strip() == "True"
    output = result.stdout + result.stderr
    assert "Server configuration" in output
    for value in ("marker-access-key-3f9a", "marker-secret-key-7c21", "marker-hcaptcha-5b8e", "marker-db-password-9d4c"):
        assert value not in output
