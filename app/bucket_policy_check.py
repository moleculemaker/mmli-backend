"""Check that no MinIO bucket grants anonymous access beyond an allowlist.

Run daily by the bucket-policy-check CronJob (chart/templates/bucket-policy-check.yaml).
Exits 0 when every bucket matches, 1 when any bucket differs from the allowlist, and 2
when the policies could not be read. Differences are logged and, when recipients are
configured, emailed.

A bucket missing from the allowlist must have no bucket policy at all, whoever its
principal names. A bucket on the allowlist must have exactly the listed policy, so a
narrowing is reported too: the allowlist has to stay a true record of what is public.

Matching is exact after sorting lists. If a MinIO upgrade starts returning an equivalent
policy in a different form (say "Principal": "*" instead of {"AWS": ["*"]}), every
allowlisted bucket is reported at once; that, rather than a real change, is the likely
cause of such an email. Update the allowlist to the new form.

Why this exists: MinIO creates buckets private, but its console and `mc anonymous set`
make a bucket public in one step. ez-specificity was found on 2026-11-05 with the canned
`public` policy, which let anyone list, read, overwrite and delete every object.

Credentials come from MINIO_ACCESS_KEY / MINIO_SECRET_KEY and need only
s3:ListAllMyBuckets, s3:GetBucketLocation and s3:GetBucketPolicy. minio-py looks up each
bucket's region before any bucket call, so GetBucketPolicy alone is denied.
"""
import json
import sys

from minio import Minio
from minio.error import S3Error

from config import MINIO_ACCESS_KEY, MINIO_SECRET_KEY, MINIO_SERVER, app_config, get_logger

log = get_logger('bucket_policy_check')


def normalize(policy):
    """Canonical JSON for a policy, or None if it grants nothing.

    MinIO returns the same policy with statements and action lists in varying order
    (staging and prod disagree), so lists are sorted before comparing.
    """
    if policy is None:
        return None
    if isinstance(policy, str):
        policy = json.loads(policy)
    if not policy.get('Statement'):
        return None

    def canon(value):
        if isinstance(value, dict):
            return {k: canon(v) for k, v in value.items()}
        if isinstance(value, list):
            return sorted((canon(v) for v in value), key=lambda v: json.dumps(v, sort_keys=True))
        return value

    return json.dumps(canon(policy), sort_keys=True)


def find_differences(actual, allowed):
    """Describe every bucket whose policy is not the allowed one.

    actual:  {bucket: policy JSON string or None} for every bucket on the server.
    allowed: {bucket: policy dict} from config.
    """
    differences = []
    for bucket in sorted(actual):
        have = normalize(actual[bucket])
        want = normalize(allowed.get(bucket))
        if have == want:
            continue
        if want is None:
            differences.append(f'{bucket}: has a bucket policy but is not on the allowlist:\n  {have}')
        elif have is None:
            differences.append(f'{bucket}: allowlisted policy is no longer set; update the allowlist')
        else:
            differences.append(f'{bucket}: policy differs from the allowlist:\n  have {have}\n  want {want}')
    for bucket in sorted(set(allowed) - set(actual)):
        differences.append(f'{bucket}: on the allowlist but no such bucket exists')
    return differences


def read_policies(client):
    policies = {}
    for bucket in client.list_buckets():
        try:
            policies[bucket.name] = client.get_bucket_policy(bucket.name)
        except S3Error as err:
            if err.code != 'NoSuchBucketPolicy':
                raise
            policies[bucket.name] = None
    return policies


def notify(recipients, server, differences):
    # Imported here so the policy logic stays importable without the email settings.
    from services.email_service import EmailService

    body = (f'The daily MinIO bucket policy check on {server} found buckets whose '
            f'policy differs from the allowlist (config.bucketPolicyCheck).\n\n'
            + '\n\n'.join(differences)
            + '\n\nIf a change was intended, update the allowlist in the chart values. '
              'Otherwise revert it, e.g. `mc anonymous set none <alias>/<bucket>`.')
    for recipient in recipients:
        try:
            EmailService().send_email(recipient, f'MinIO bucket policy check failed on {server}', body)
        except Exception as ex:
            log.error(f'Failed to email {recipient}: {ex}')


def main():
    settings = app_config['bucketPolicyCheck']
    allowed = settings.get('allowedAnonymous') or {}
    recipients = settings.get('recipients') or []

    client = Minio(MINIO_SERVER, access_key=MINIO_ACCESS_KEY, secret_key=MINIO_SECRET_KEY, secure=False)
    # The comparison is inside the try too: a policy that is not valid JSON fails in
    # normalize(), and that must take the same emailed exit-2 path as a failed read.
    try:
        actual = read_policies(client)
        differences = find_differences(actual, allowed)
    except Exception as ex:
        log.error(f'Could not read or parse bucket policies from {MINIO_SERVER}: {ex}')
        notify(recipients, MINIO_SERVER, [f'The check could not read or parse bucket policies: {ex}'])
        return 2

    if not differences:
        log.info(f'All {len(actual)} buckets match the allowlist')
        return 0
    for difference in differences:
        log.error(difference)
    notify(recipients, MINIO_SERVER, differences)
    return 1


if __name__ == '__main__':
    sys.exit(main())
