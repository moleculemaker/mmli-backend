"""A CLEAN submission with no sequences must be rejected, not run.

The legacy route skips schema validation, and the CLEAN command builder only checked the
upper bound (20 sequences). Job 8e8bca70 (2026-07-07) was submitted with
`{"input_fasta": {}}`: it ran CLEAN on an empty FASTA, which crashed, and was marked
completed with no results.
"""
import json

import pytest

SEQUENCE = {"header": "seq1", "sequence": "MKTAYIAKQR", "DNA_sequence": ""}


@pytest.mark.parametrize("input_fasta", [{}, []], ids=["empty-object", "empty-list"])
def test_legacy_rejects_no_sequences(client, input_fasta):
    resp = client.post('/clean/jobs', json={'job_info': json.dumps({'input_fasta': input_fasta})})

    assert resp.status_code == 400
    assert 'at least one' in resp.json()['detail']
    assert client.k8s_jobs == []


def test_legacy_accepts_one_sequence(client):
    resp = client.post('/clean/jobs', json={'job_info': json.dumps({'input_fasta': [SEQUENCE]})})

    assert resp.status_code == 201, resp.text
    assert len(client.k8s_jobs) == 1
