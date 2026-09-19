"""Tests for the SimpleFold sequence-length bound.

Why this exists: `ml-simplefold`'s VRAM grows with chain length until it exhausts the
shared GPU. At 1022 residues it OOMs running ALONE -- and 1022 is exactly what the MEP
form used to accept, because that is ESM-2's context limit and has nothing to do with
folding. Those submissions fail in production today.

The bound is enforced server-side rather than only in the browser because the browser is
not the only client. Both `/{job_type}/jobs` and `/v1/...` reach `job_builder.prepare_job`
before `kubejob_service.create_job` (two callers of each, checked), and `coordinator.py`
creates subjobs by POSTing to the same router, so the check in `prepare_job` is the one
choke point every path crosses. The tests below drive `prepare_job` and both routers, not
a helper, so the wiring (config key, job type, check-before-upload order, error shape per
surface) is what is asserted.
"""
import json
import pathlib

import pytest
import yaml
from fastapi import HTTPException

from models.enums import JobType
from routers.v1.problems import INVALID_INPUT
from services import job_builder
from services.fasta import count_fasta_residues

# The limit the checked-in config carries. Tests that need a specific number read it
# from the config rather than restating it, so raising the limit does not go red here.
SIMPLEFOLD = JobType.ML_SIMPLEFOLD
LIMIT = job_builder.app_config['kubernetes_jobs'][SIMPLEFOLD]['maxResidues']


def _fasta(n: int) -> str:
    return '>h\n' + 'A' * n


def _prepare(fake_minio, fasta, job_type=SIMPLEFOLD):
    return job_builder.prepare_job(
        job_type=job_type, job_id='j1', job_info=json.dumps({'fasta': fasta}),
        service=fake_minio,
    )


class TestCountFastaResidues:
    @pytest.mark.parametrize("fasta,expected,why", [
        (">h\nABCDE", 5, "the shape the frontend sends"),
        (">h\nABC\nDE", 5, "pasted FASTA is line-wrapped"),
        ("ABCDE", 5, "bare sequence, no header"),
        (">h\nABCDE*", 5, "trailing stop codon is not a residue"),
        (">h\nABC\n\nDE\n", 5, "blank lines"),
        (">h\n ABC DE ", 5, "internal and surrounding whitespace"),
        ("", 0, "empty input must not raise"),
        (">h\n", 0, "header with no body"),
    ])
    def test_single_record(self, fasta, expected, why):
        assert count_fasta_residues(fasta) == expected, why

    @pytest.mark.parametrize("fasta", [">a\nABC\n>b\nABCDEFG", ">b\nABCDEFG\n>a\nABC"])
    def test_multi_record_takes_the_longest_not_the_total(self, fasta):
        """Peak GPU memory is driven by the largest single chain, so summing would reject
        many short chains that fold fine, and averaging would admit one that cannot."""
        assert count_fasta_residues(fasta) == 7

    @pytest.mark.parametrize("fasta,expected,why", [
        (">h\nAB--CDE", 7, "alignment gaps reach the model as positions"),
        (">h\n1 ABCDE 6", 7, "position digits reach the model as positions"),
        (">h\n;comment\nABCDE", 13, "a ';' line after the header is body to Biopython"),
    ])
    def test_non_letter_body_characters_count(self, fasta, expected, why):
        """SimpleFold reads the file with Biopython and hands `str(record.seq)` to the
        model unchanged, so every non-whitespace body character costs GPU memory. This
        pins that the counter measures what the GPU sees, not what a residue alphabet
        would admit -- counting letters only would under-count."""
        assert count_fasta_residues(fasta) == expected, why


class TestPrepareJobEnforcesTheBound:
    def test_under_the_bound_uploads_the_fasta(self, fake_minio):
        _prepare(fake_minio, _fasta(LIMIT - 1))
        assert fake_minio.get_file(SIMPLEFOLD, '/j1/in/input.fasta') is not None

    def test_exactly_at_the_bound_passes(self, fake_minio):
        """The limit is inclusive. LIMIT is a length we intend to accept, so an
        off-by-one here would silently narrow the product by one residue."""
        _prepare(fake_minio, _fasta(LIMIT))

    def test_over_the_bound_is_422_and_nothing_is_uploaded(self, fake_minio):
        with pytest.raises(job_builder.InputTooLarge) as e:
            _prepare(fake_minio, _fasta(LIMIT + 1))
        # 422, not 400: the request is well-formed, the input is simply too large for
        # the hardware. It is also not transient, so a client must not retry it.
        assert e.value.status_code == 422
        assert e.value.pointer == '/fasta'
        # The check runs before MinIO, so a rejected submission leaves nothing behind.
        assert fake_minio.objects == {}

    def test_the_error_names_both_numbers(self, fake_minio):
        """A bare rejection sends the user back to guess. The message has to say what
        they submitted and what the limit is, because the limit is not in the UI copy
        on the API path."""
        with pytest.raises(HTTPException) as e:
            _prepare(fake_minio, _fasta(1022))
        assert '1022' in str(e.value.detail) and str(LIMIT) in str(e.value.detail)

    def test_the_1022_case_that_fails_in_production_is_rejected(self, fake_minio):
        """The regression this bound exists for: 1022 residues is what the MEP form
        allowed, and simplefold OOMs on it even with the GPU to itself."""
        with pytest.raises(HTTPException):
            _prepare(fake_minio, _fasta(1022))

    @pytest.mark.parametrize("bad", [None, ['>h', 'ABC'], 7])
    def test_a_non_string_fasta_is_400_not_500(self, fake_minio, bad):
        """The legacy router has no schema validation, so these reach prepare_job as-is
        and used to AttributeError on `.encode`."""
        with pytest.raises(HTTPException) as e:
            _prepare(fake_minio, bad)
        assert e.value.status_code == 400

    def test_a_missing_config_key_fails_loudly(self, fake_minio, monkeypatch):
        """No fallback constant: an operator who misspells `maxResidues` in the chart
        must find out on the first submission, not deploy a limit that never applies."""
        config = json.loads(json.dumps(job_builder.app_config))  # deep copy
        del config['kubernetes_jobs'][SIMPLEFOLD]['maxResidues']
        monkeypatch.setattr(job_builder, 'app_config', config)
        with pytest.raises(KeyError):
            _prepare(fake_minio, _fasta(1))

    def test_the_mep_job_is_not_bounded(self, fake_minio):
        """770 residues (APP770) is the longest sequence in 173 historical MEP
        submissions. Only the structure half is bounded; the MEP job itself must still
        accept anything up to ESM-2's 1022, so this drives the MEP branch of
        prepare_job, not the counter."""
        prepared = job_builder.prepare_job(
            job_type=JobType.CLEANDB_MEPESM, job_id='m1',
            job_info=json.dumps({'sequence': 'A' * 770}), service=fake_minio,
        )
        assert {'name': 'CLEANDB_INPUT_SEQUENCE', 'value': 'A' * 770} in prepared.environment


# The pod does not read `app/cfg/config.yaml` as checked in: the chart renders a
# ConfigMap from `{{ .Values.config }}` and mounts it over that path, with each
# environment's values file coalesced onto `chart/values.yaml`. The hard index above
# turns a missing key into a 500 on the first submission, so every deployed
# combination is checked for the key here.
_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_BASE_VALUES = 'chart/values.yaml'
DEPLOYED_VALUES_FILES = [
    _BASE_VALUES,
    'chart/values.prod.yaml',
    'chart/values.staging.yaml',
    'chart/values.mmli2.prod.yaml',
    'chart/values.mmli2.staging.yaml',
]


def _coalesce(base, overlay):
    """Merge the way Helm does: maps merge key by key, anything else the overlay wins."""
    merged = dict(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _coalesce(merged[key], value)
        else:
            merged[key] = value
    return merged


def _deployed_config(values_file):
    base = yaml.safe_load((_REPO_ROOT / _BASE_VALUES).read_text())
    if values_file == _BASE_VALUES:
        return base['config']
    overlay = yaml.safe_load((_REPO_ROOT / values_file).read_text())
    return _coalesce(base, overlay)['config']


class TestEveryDeployedConfigCarriesTheBound:
    @pytest.mark.parametrize("values_file", DEPLOYED_VALUES_FILES)
    def test_the_check_resolves_its_key_in_every_environment(
            self, values_file, fake_minio, monkeypatch):
        config = _deployed_config(values_file)
        monkeypatch.setattr(job_builder, 'app_config', config)
        limit = config['kubernetes_jobs'][SIMPLEFOLD]['maxResidues']
        _prepare(fake_minio, _fasta(limit))
        with pytest.raises(job_builder.InputTooLarge):
            _prepare(fake_minio, _fasta(limit + 1))


class TestEachSurfaceReportsItInItsOwnShape:
    def test_legacy_returns_422_with_a_detail(self, client):
        resp = client.post(f'/{SIMPLEFOLD}/jobs',
                           json={'job_info': json.dumps({'fasta': _fasta(LIMIT + 1)})})
        assert resp.status_code == 422
        assert str(LIMIT) in resp.json()['detail']
        assert client.k8s_jobs == []

    def test_v1_returns_an_invalid_input_problem_pointing_at_fasta(self, client):
        """The MCP adapter relays this body verbatim so an agent can correct its call.
        A bare HTTPException would render as `about:blank` with no pointer, which the
        adapter cannot classify; this pins the same shape schema failures use."""
        resp = client.post(f'/v1/tools/{SIMPLEFOLD}/jobs', json={'fasta': _fasta(LIMIT + 1)})
        assert resp.status_code == 422
        assert resp.headers['content-type'].startswith('application/problem+json')
        body = resp.json()
        assert body['type'] == INVALID_INPUT
        assert body['errors'][0]['pointer'] == '/fasta'
        assert str(LIMIT) in body['detail']
        assert client.k8s_jobs == []

    def test_v1_accepts_the_bound_itself(self, client):
        resp = client.post(f'/v1/tools/{SIMPLEFOLD}/jobs', json={'fasta': _fasta(LIMIT)})
        assert resp.status_code == 201
        assert len(client.k8s_jobs) == 1

    def test_the_descriptor_advertises_the_bound(self, client):
        """One source of truth for both the frontend's own check and API clients: the
        limit a client can read before submitting is the limit the check enforces."""
        assert client.get(f'/v1/tools/{SIMPLEFOLD}').json()['x-max-residues'] == LIMIT
        # Unbounded tools publish null rather than omitting the key.
        assert client.get('/v1/tools/somn').json()['x-max-residues'] is None
