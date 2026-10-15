"""Tests pinning which SimpleFold model variant actually runs.

The defect these exist for: the image tag and the model the backend asks for are two
separate declarations, and the Python one silently outranks the YAML one. `create_job`
passes a command, which overrides the container's ENTRYPOINT, so the image's own
`predict-simplefold.sh` (which selects simplefold_1.6B) never executes.

That is how `cee8751` -- a commit titled "chore: alphabetize kubernetes_job in
values.yaml" -- bumped the image from `version1` to `version1.6B` without changing what
ran. Every job has used simplefold_100M since the tool was added in `6903e4a`, and
nothing anywhere said so.

Nothing here asserts that 100M is the RIGHT variant; that is a science call. These only
ensure the choice is visible, single-sourced, and cannot change by accident.
"""
import copy
import json

import pytest
import yaml
from fastapi import HTTPException

from models.enums import JobType
from services import job_builder
from services.job_builder import SIMPLEFOLD_MODELS

from conftest import BASE_VALUES, DEPLOYED_VALUES_FILES, REPO_ROOT, deployed_config

# What local runs and this test suite load. The chart replaces it entirely when deployed,
# so it is not one of DEPLOYED_VALUES_FILES, but it declares the same job and has to be
# checked the same way.
IN_REPO_CONFIG = "app/cfg/config.yaml"


def _simplefold_config(source):
    if source == IN_REPO_CONFIG:
        doc = yaml.safe_load((REPO_ROOT / source).read_text())
    else:
        doc = deployed_config(source)
    return doc["kubernetes_jobs"]["ml-simplefold"]


def _prepare(monkeypatch, fake_minio, **overrides):
    """Run prepare_job for SimpleFold against the loaded config with the ml-simplefold
    block's `overrides` applied. A value of None removes the key."""
    cfg = copy.deepcopy(job_builder.app_config)
    block = cfg["kubernetes_jobs"]["ml-simplefold"]
    for key, value in overrides.items():
        if value is None:
            block.pop(key, None)
        else:
            block[key] = value
    monkeypatch.setattr(job_builder, "app_config", cfg)
    return job_builder.prepare_job(
        JobType.ML_SIMPLEFOLD, "job-1", json.dumps({"fasta": ">a\nMKV"}), fake_minio)


class TestTheConfiguredVariantIsWhatRuns:
    """The behavior the whole change exists for: the key in config is the model on the
    command line. This is the test that fails if the key name drifts between the YAML
    and the code, or if a variant is hardcoded into the command again."""

    def test_prepare_job_passes_the_configured_model(self, fake_minio, monkeypatch):
        prepared = _prepare(monkeypatch, fake_minio, simplefoldModel="simplefold_1.6B")
        assert " --simplefold_model simplefold_1.6B " in prepared.command

    def test_a_missing_key_fails_instead_of_picking_a_variant(self, fake_minio, monkeypatch):
        """A silent default is how the variant became invisible in the first place."""
        with pytest.raises(KeyError):
            _prepare(monkeypatch, fake_minio, simplefoldModel=None)
        assert fake_minio.objects == {}

    def test_a_variant_the_image_does_not_ship_fails_before_launch(self, fake_minio, monkeypatch):
        """Otherwise the job is scheduled onto the shared GPU and dies minutes later with
        FileNotFoundError on artifacts/simplefold_700M.ckpt."""
        with pytest.raises(HTTPException) as raised:
            _prepare(monkeypatch, fake_minio, simplefoldModel="simplefold_700M")
        assert "simplefold_700M" in raised.value.detail
        # The config is checked before the FASTA is uploaded, so a submission rejected
        # for a misconfigured variant leaves nothing behind in MinIO.
        assert fake_minio.objects == {}


class TestEveryDeploymentDeclaresAShippedVariant:
    @pytest.mark.parametrize("source", DEPLOYED_VALUES_FILES + [IN_REPO_CONFIG])
    def test_the_model_is_declared_and_exists_in_the_image(self, source):
        """Absent or not shipped in the image, job creation fails (above), so every
        SimpleFold job in that cluster would fail. Caught here, per cluster, before a
        deploy rather than after."""
        assert _simplefold_config(source).get("simplefoldModel") in SIMPLEFOLD_MODELS

    def test_local_runs_use_the_same_model_as_the_default_deployment(self):
        """Overlays may legitimately differ per cluster, but app/cfg/config.yaml and the
        base chart are both 'the default', and a split between them means local runs
        and deployed runs use different models -- the same class of split this change
        exists to close."""
        assert (_simplefold_config(IN_REPO_CONFIG)["simplefoldModel"]
                == _simplefold_config(BASE_VALUES)["simplefoldModel"])
