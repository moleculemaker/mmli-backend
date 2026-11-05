"""The daily MinIO bucket policy check must flag any anonymous access not on the allowlist.

The policies below are the real ones read from the mmli2 servers on 2026-11-05.
"""
import json
from types import SimpleNamespace

import pytest
from minio.error import S3Error

import bucket_policy_check as check

# Canned `public` policy as staging returned it ...
PUBLIC_STAGING = ('{"Statement":[{"Action":["s3:ListBucketMultipartUploads","s3:GetBucketLocation","s3:ListBucket"],'
                  '"Effect":"Allow","Principal":{"AWS":["*"]},"Resource":["arn:aws:s3:::ez-specificity"]},'
                  '{"Action":["s3:PutObject","s3:AbortMultipartUpload","s3:DeleteObject","s3:GetObject","s3:ListMultipartUploadParts"],'
                  '"Effect":"Allow","Principal":{"AWS":["*"]},"Resource":["arn:aws:s3:::ez-specificity/*"]}],"Version":"2012-10-17"}')
# ... and as prod returned it: same policy, different action order.
PUBLIC_PROD = ('{"Statement":[{"Action":["s3:ListBucket","s3:ListBucketMultipartUploads","s3:GetBucketLocation"],'
               '"Effect":"Allow","Principal":{"AWS":["*"]},"Resource":["arn:aws:s3:::ez-specificity"]},'
               '{"Action":["s3:AbortMultipartUpload","s3:DeleteObject","s3:GetObject","s3:ListMultipartUploadParts","s3:PutObject"],'
               '"Effect":"Allow","Principal":{"AWS":["*"]},"Resource":["arn:aws:s3:::ez-specificity/*"]}],"Version":"2012-10-17"}')

COMPLEXES_ONLY = {
    "Version": "2012-10-17",
    "Statement": [{"Effect": "Allow", "Principal": {"AWS": ["*"]}, "Action": ["s3:GetObject"],
                   "Resource": ["arn:aws:s3:::ez-specificity/*/out/*/out/complexes/*"]}],
}


def _s3_error(code):
    return S3Error(code, 'msg', 'resource', 'request-id', 'host-id', None)


class FakeClient:
    def __init__(self, policies):
        # bucket -> policy JSON string, None (no policy), or an exception to raise
        self.policies = policies

    def list_buckets(self):
        return [SimpleNamespace(name=name) for name in self.policies]

    def get_bucket_policy(self, bucket):
        policy = self.policies[bucket]
        if isinstance(policy, Exception):
            raise policy
        if policy is None:
            raise _s3_error('NoSuchBucketPolicy')
        return policy


def test_the_same_policy_in_a_different_order_normalizes_equal():
    assert check.normalize(PUBLIC_STAGING) == check.normalize(PUBLIC_PROD)


def test_an_empty_statement_list_counts_as_no_policy():
    # staging's aceretro bucket returns exactly this.
    assert check.normalize('{"Statement":[],"Version":"2012-10-17"}') is None


def test_the_public_policy_found_on_ez_specificity_is_flagged():
    differences = check.find_differences({'ez-specificity': PUBLIC_PROD}, {'ez-specificity': COMPLEXES_ONLY})
    assert len(differences) == 1
    assert 'differs from the allowlist' in differences[0]


def test_any_policy_on_a_bucket_not_on_the_allowlist_is_flagged():
    differences = check.find_differences({'clean': json.dumps(COMPLEXES_ONLY)}, {})
    assert differences and 'not on the allowlist' in differences[0]


def test_private_buckets_not_on_the_allowlist_pass():
    assert check.find_differences({'clean': None, 'aceretro': '{"Statement":[]}'}, {}) == []


def test_an_allowlisted_policy_that_was_removed_is_reported():
    differences = check.find_differences({'ez-specificity': None}, {'ez-specificity': COMPLEXES_ONLY})
    assert differences and 'no longer set' in differences[0]


def test_an_allowlisted_bucket_that_does_not_exist_is_reported():
    differences = check.find_differences({}, {'ez-specificity': COMPLEXES_ONLY})
    assert differences and 'no such bucket' in differences[0]


def test_no_such_bucket_policy_is_read_as_private():
    assert check.read_policies(FakeClient({'clean': None})) == {'clean': None}


def test_any_other_read_error_is_raised_not_treated_as_private():
    with pytest.raises(S3Error):
        check.read_policies(FakeClient({'clean': _s3_error('AccessDenied')}))


@pytest.fixture
def run(monkeypatch):
    """Run main() against a fake server; returns (exit code, emailed differences)."""
    def _run(policies, allowed, recipients=('ops@example.org',)):
        sent = []
        monkeypatch.setitem(check.app_config, 'bucketPolicyCheck',
                            {'allowedAnonymous': allowed, 'recipients': list(recipients)})
        monkeypatch.setattr(check, 'Minio', lambda *a, **k: FakeClient(policies))
        monkeypatch.setattr(check, 'notify', lambda recipients, server, differences: sent.append(differences))
        return check.main(), sent
    return _run


def test_main_exits_0_and_sends_nothing_when_everything_matches(run):
    code, sent = run({'clean': None, 'ez-specificity': json.dumps(COMPLEXES_ONLY)}, {'ez-specificity': COMPLEXES_ONLY})
    assert (code, sent) == (0, [])


def test_main_exits_1_and_notifies_on_a_difference(run):
    code, sent = run({'ez-specificity': PUBLIC_STAGING}, {'ez-specificity': COMPLEXES_ONLY})
    assert code == 1
    assert len(sent) == 1 and 'ez-specificity' in sent[0][0]


def test_main_exits_2_and_notifies_when_policies_cannot_be_read(run):
    code, sent = run({'clean': _s3_error('AccessDenied')}, {})
    assert code == 2
    assert len(sent) == 1 and 'could not read' in sent[0][0]


def test_notify_sends_one_email_per_recipient_and_survives_a_failure(monkeypatch):
    from services import email_service

    calls = []

    def fake_send(self, recipients=None, email_subject='', message_body='', template_file=None):
        calls.append(recipients)
        if recipients == 'a@example.org':
            raise OSError('mail server down')

    monkeypatch.setattr(email_service.EmailService, 'send_email', fake_send)
    check.notify(['a@example.org', 'b@example.org'], 'minio:9000', ['x'])
    assert calls == ['a@example.org', 'b@example.org']
