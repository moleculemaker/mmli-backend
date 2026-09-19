"""Tests for KubeEventWatcher.send_notification_email.

The bug this file starts from: `cleandb-mepesm` (the Mutation Effect Prediction tool,
MEP-ESM) matched none of the per-tool branches, so it fell through to the catch-all

    elif job_type in JobTypes:   # "OED & CLEANDB jobs are very fast"
        return

and no completion email was ever sent -- while the frontend's submission form offers the
user an address field and a checkbox reading "Agree to receive email notifications".
The rationale in that comment was also simply untrue by the time it was read: a MEP-ESM
run waits several minutes for the single GPU node before it computes.

The failure mode is quiet by construction. Nothing raises, nothing retries, and the only
trace is a `log.warning` on a line that looks deliberate, so the shape worth pinning is
not "does an email go out" in general but "does THIS job type reach the send". The
catch-all makes every future tool fail the same way by default, which is why the enum is
partitioned explicitly below into the types that notify and the types deliberately kept
silent -- a fix that widens the dispatch far enough to notify ezspec's intermediary
subjobs would be its own bug, and a tool added to neither list is the original bug again.
"""
import pathlib

import pytest
import yaml

from models.enums import JobStatus, JobType
from services import kubejob_service
from services.kubejob_service import KubeEventWatcher


class _Job:
    """Minimal stand-in for the Job row the watcher hands to the notifier."""

    def __init__(self, job_id="job-abc123", email="scientist@illinois.edu"):
        self.job_id = job_id
        self.email = email


class _FakeEmailService:
    def __init__(self):
        self.sent = []

    def send_email(self, recipients, subject, body):
        self.sent.append({"to": recipients, "subject": subject, "body": body})


def _watcher():
    """Build a watcher without its __init__ side effects.

    conftest has already replaced KubeEventWatcher.__init__ with a no-op (it otherwise
    opens a DB connection and starts an endless watch thread), so constructing one here
    is safe; only the one collaborator the notifier still touches is supplied.

    No MinIO stand-in any more: idempotency moved off the `{job_id}/email-sent` object
    and onto `job.notified_at`, which the caller owns. The notifier has no memory, so
    "already notified" is not a state this file can be in - see
    tests/test_watcher_notifications.py for the gate itself.
    """
    w = KubeEventWatcher()
    w.email_service = _FakeEmailService()
    return w


def _notify(watcher, job_type, phase=JobStatus.COMPLETED, job=None):
    """Call the notifier the way the watch loop does.

    `job_type` is passed as a plain `str` on purpose: the real caller reads it out of the
    Kubernetes label `jobType`, never as a JobType member. JobType subclasses str, so the
    equality checks work either way -- but only the string reproduces the live path.
    """
    job = job or _Job()
    watcher.send_notification_email(job.job_id, str(job_type), job, phase)
    return watcher.email_service.sent


class TestMepEsmIsNotified:
    def test_completion_sends_an_email(self):
        """THE REGRESSION. Before the fix this list was empty: cleandb-mepesm reached the
        catch-all and returned without sending."""
        w = _watcher()
        sent = _notify(w, JobType.CLEANDB_MEPESM)
        assert len(sent) == 1

    def test_completion_email_links_to_the_result_page(self):
        """The route is 'effect-prediction/result/:id' -- singular 'result', unlike the
        '/results/' every neighbouring branch uses. Getting it wrong sends a live link to
        a 404, which is worse than the silence this replaces, so the shape is pinned."""
        w = _watcher()
        sent = _notify(w, JobType.CLEANDB_MEPESM)
        assert f"/effect-prediction/result/{_Job().job_id}" in sent[0]["body"]

    def test_failure_sends_an_email(self):
        w = _watcher()
        sent = _notify(w, JobType.CLEANDB_MEPESM, phase=JobStatus.ERROR)
        assert len(sent) == 1
        assert "failed" in sent[0]["subject"].lower()

    def test_a_send_reports_itself_so_the_caller_can_record_it(self):
        """The watcher re-reconciles every listed Job on each pass, so a completed job
        reaches the notifier repeatedly. What stops a re-send is the caller writing
        job.notified_at, and it only does that when this returns True."""
        w = _watcher()
        assert w.send_notification_email("job-abc123", str(JobType.CLEANDB_MEPESM),
                                         _Job(), JobStatus.COMPLETED) is True

    def test_a_job_with_no_address_sends_nothing(self):
        w = _watcher()
        assert _notify(w, JobType.CLEANDB_MEPESM, job=_Job(email=None)) == []

    def test_a_job_with_no_address_reports_nothing_sent(self):
        """False rather than True matters here: True would have the caller stamp
        notified_at on a job that was never mailed."""
        w = _watcher()
        assert w.send_notification_email("job-abc123", str(JobType.CLEANDB_MEPESM),
                                         _Job(email=None), JobStatus.COMPLETED) is False

    def test_a_non_terminal_phase_reports_nothing_sent(self):
        """The clause the caller leans on most often, and the only one whose failure is
        silent: a truthy return for 'processing' would stamp notified_at on a running
        job and lose its completion email. The watcher tests cover this through a stub,
        which cannot catch the real dispatch drifting."""
        w = _watcher()
        assert w.send_notification_email("job-abc123", str(JobType.CLEANDB_MEPESM),
                                         _Job(), JobStatus.PROCESSING) is False

    def test_a_failed_send_reports_nothing_sent_so_it_stays_retryable(self):
        def _boom(*args, **kwargs):
            raise RuntimeError("SMTP unavailable")

        w = _watcher()
        w.email_service.send_email = _boom

        assert w.send_notification_email("job-abc123", str(JobType.CLEANDB_MEPESM),
                                         _Job(), JobStatus.COMPLETED) is False


# Every JobType must be listed in exactly one of these two sets.
#
# The dispatch's catch-all (`elif job_type in JobTypes`) returns silently, and `JobTypes`
# is derived from the enum -- `[str(job_type) for job_type in JobType]` in
# `app/models/enums.py` -- so every member satisfies it by construction. A tool added to
# the enum and forgotten in the dispatch therefore never raises and never emails; it just
# logs a warning that reads as intentional. That is precisely how CLEANDB_MEPESM went
# unnotified, and it is why "no known type raises" cannot be the safety net: the
# `raise ValueError` below the catch-all is unreachable for any JobType member.
#
# The net is instead this explicit partition. Adding a JobType without classifying it
# here fails `test_every_job_type_is_classified`, which forces the choice to be made
# deliberately rather than defaulted into silence.

NOTIFYING_JOB_TYPES = [
    JobType.ACERETRO,
    JobType.CLEAN,
    JobType.CLEANDB_MEPESM,
    JobType.CRISPR_COPIES,
    JobType.EZ_SPECIFICITY,
    JobType.MOLLI,
    JobType.MUTAGENESIS,
    JobType.NOVOSTOIC_OPTSTOIC,
    JobType.NOVOSTOIC_PATHWAYS,
    JobType.NOVOSTOIC_ENZRANK,
    JobType.NOVOSTOIC_DGPREDICTOR,
    JobType.OED_CHEMINFO,
    JobType.REACTIONMINER,
    JobType.SOMN,
]

DELIBERATELY_SILENT_JOB_TYPES = [
    JobType.EZSPEC_UNIDOCK,     # intermediary step of an ez-specificity run
    JobType.EZSPEC_INFERENCE,   # ditto; the parent job is what the user waits on
    JobType.OED_DLKCAT,
    JobType.OED_UNIKP,
    JobType.OED_CATPRED,
    # ML_SIMPLEFOLD is silent because it is the COMPANION half of a cleandb-mepesm
    # submission, not because it lacks a frontend (the dispatch's own comment says
    # "no frontend URL yet"; that is stale -- the MEP result page renders its
    # structure). CLEANDB-frontend's effect-prediction submit handler creates the
    # simplefold job and then the MEP job, passing the SAME email to both, so giving
    # this one a branch would send the user two emails for one submission. The
    # MEP-ESM notification covers the pair.
    JobType.ML_SIMPLEFOLD,
    JobType.DEFAULT,            # example jobs
    # CHEMSCRAPER is listed here to describe today's behavior, but it is NOT a deliberate
    # skip -- it is the same live defect this PR fixes for MEP-ESM, and it needs its own
    # change. Do not read its presence in this list as the question being settled.
    #
    # Chemscraper does have an emailing path (chemscraper_service.py:291/:302, success
    # and failure), but it hangs off `POST /chemscraper/analyze`, which is marked
    # `deprecated=True` at routers/chemscraper.py:31 and which the frontend never calls.
    # What the frontend actually calls is the generic route:
    #
    #   configuration.component.html:105   email input -> `userEmail`
    #   configuration.component.ts:114     userEmail -> requestBody.user_email
    #   configuration.component.ts:90      -> chemscraper.service.ts:44 analyzeDocument()
    #   chemscraper.service.ts:47          createJobJobTypeJobsPost('chemscraper', {email})
    #                                      == POST /chemscraper/jobs
    #   routers/job.py:62 -> job_builder.py:134 -> a real Kubernetes Job
    #
    # The watcher observes that Job, reaches this dispatch, finds no CHEMSCRAPER branch,
    # and returns. So a user who typed an address into the chemscraper form is never
    # notified -- exactly the MEP-ESM bug, still open. (Verified against a fresh clone of
    # moleculemaker/chemscraper-frontend @ 697af1a; the method is named analyzeDocument
    # after the endpoint it no longer uses, which is what hid this.)
    JobType.CHEMSCRAPER,
]


class TestEveryJobTypeIsClassified:
    """The safety net proper: a new tool cannot default into silence unnoticed."""

    def test_every_job_type_is_classified(self):
        classified = set(NOTIFYING_JOB_TYPES) | set(DELIBERATELY_SILENT_JOB_TYPES)
        unclassified = set(JobType) - classified
        assert not unclassified, (
            f"JobType(s) {sorted(str(j) for j in unclassified)} are in neither list. The "
            f"dispatch will skip them silently -- add a branch in send_notification_email "
            f"and list them under NOTIFYING_JOB_TYPES, or record the skip as deliberate."
        )

    def test_the_two_sets_are_disjoint(self):
        assert not set(NOTIFYING_JOB_TYPES) & set(DELIBERATELY_SILENT_JOB_TYPES)


class TestNotifyingTypesSendAndResolveTheirConfigKey:
    """Doubles as config-key coverage: each branch reads `app_config['<tool>_frontend_url']`
    unconditionally, so a key missing from `app/cfg/config.yaml` raises `KeyError` here.
    This is what caught the absent `reactionminer_frontend_url`.

    Note the file: conftest points CONFIG_FILEPATH at `app/cfg/config.yaml`, which is NOT
    what any cluster reads. The chart values are covered separately below."""

    @pytest.mark.parametrize("job_type", NOTIFYING_JOB_TYPES, ids=str)
    def test_one_email_is_sent(self, job_type):
        assert len(_notify(_watcher(), job_type)) == 1


# The class above reads `app/cfg/config.yaml`. In every deployed environment that file is
# shadowed: `chart/templates/deployment.yaml` mounts the ConfigMap over
# `/code/app/cfg/config.yaml` with `subPath: config.yaml`, and
# `chart/templates/configmap.yaml` renders that ConfigMap from `{{ .Values.config }}`.
# The keys the dispatch reads therefore have to exist in the CHART values, and nothing
# checked those -- which is why adding `cleandb_frontend_url` here meant editing eight
# files by hand with no test to say whether one had been missed.
#
# What a pod sees is base + overlay, because Helm coalesces an overlay onto
# `chart/values.yaml` map key by map key. That merge is reconstructed below, so a key
# absent from BOTH raises the same `KeyError` inside `send_notification_email` that the
# pod would raise -- where `run()`'s `except Exception` swallows it and that one tool
# silently stops notifying in that one cluster.
#
# A key present in the base but missing from an overlay is NOT this failure: the base
# value is inherited, so the pod gets a working link to the wrong environment. That is a
# real problem (`ezspecificity_frontend_url` is absent from both mmli1 files today, so
# mmli1 EZspecificity links at mmli2 staging) but it is a wrong-value bug rather than a
# missing-key one, and it is not what this test is for.

_REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
_BASE_VALUES = "chart/values.yaml"

# `chart/values.local.yaml` is deliberately not here. It keeps its `*_frontend_url` keys
# under `controller:` rather than `config:`, where the ConfigMap template cannot reach
# them (see the `TODO: fit this into new config structure` in that file), so merging it
# would only re-assert the base. Restructuring it is its own change.
DEPLOYED_VALUES_FILES = [
    _BASE_VALUES,                       # the fallback every overlay inherits from
    "chart/values.prod.yaml",           # mmli1
    "chart/values.staging.yaml",        # mmli1
    "chart/values.mmli2.prod.yaml",
    "chart/values.mmli2.staging.yaml",
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
    """The `config:` block a pod deployed with `values_file` would actually be handed."""
    base = yaml.safe_load((_REPO_ROOT / _BASE_VALUES).read_text())
    if values_file == _BASE_VALUES:
        return base["config"]
    overlay = yaml.safe_load((_REPO_ROOT / values_file).read_text())
    return _coalesce(base, overlay)["config"]


class TestEveryEnvironmentsConfigMapResolvesEveryKey:
    """The same coverage as above, against the config each cluster is really given."""

    @pytest.mark.parametrize("values_file", DEPLOYED_VALUES_FILES)
    def test_every_notifying_type_resolves_its_frontend_url(self, values_file, monkeypatch):
        monkeypatch.setattr(kubejob_service, "app_config", _deployed_config(values_file))

        missing = []
        for job_type in NOTIFYING_JOB_TYPES:
            try:
                _notify(_watcher(), job_type)
            except KeyError as missing_key:
                missing.append(f"{job_type} reads {missing_key}")

        assert not missing, (
            f"{values_file} (coalesced onto {_BASE_VALUES}) is missing config keys the "
            f"dispatch reads unconditionally, so send_notification_email will raise "
            f"KeyError in that environment and the email will be silently lost:\n  "
            + "\n  ".join(missing)
        )


class TestTheSilentTypesStaySilent:
    """These types send nothing today; the list above is what says which of those
    silences are deliberate and which are open (CHEMSCRAPER is open). This class asserts
    only the outcome, so do not read a passing test here as the question being settled.
    What it does pin is that widening the dispatch to reach MEP-ESM did not start
    notifying any of them.

    Note the catch-all is not the route for all of them: ML_SIMPLEFOLD and DEFAULT return
    earlier in the dispatch and never reach it."""

    @pytest.mark.parametrize("job_type", DELIBERATELY_SILENT_JOB_TYPES, ids=str)
    def test_no_email_is_sent(self, job_type):
        assert _notify(_watcher(), job_type) == []


class TestAnUnknownJobTypeStillRaises:
    def test_an_unknown_job_type_still_raises(self):
        """The guard itself is worth keeping: it is what surfaces a label the backend has
        never heard of, rather than silently dropping it. Note it is reachable only for a
        label that is not a JobType at all -- every enum member is absorbed by the
        catch-all above it."""
        with pytest.raises(ValueError):
            _notify(_watcher(), "some-tool-that-does-not-exist")
