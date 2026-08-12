"""Tests for KubeEventWatcher._derive_phase.

This is the function the stranded-job bug came down to. The watcher previously decided a
job had succeeded by testing `conditions[0].type == 'SuccessCriteriaMet'`, which fails in
three separate ways:

  * 'SuccessCriteriaMet' only exists on Kubernetes 1.31 and later. Every successful Job
    gets 'Complete', on every version, and 'Complete' was never checked at all.
  * The index is fixed, so the test depends on the order the Job controller happened to
    append conditions in.
  * A terminal event that matched neither branch fell through to an `else` that logged
    '>> Skipped job update' and moved on.

None of that is recoverable after the fact: a watch stream replays nothing older than the
resourceVersion it reconnects at, so a missed completion leaves the row at whatever phase
it last held -- 'queued', if no event for it was ever processed -- permanently.

The cases below are pinned to that history. Each one marked NOT DETECTED BEFORE is a shape
the old logic got wrong, and is the reason this file exists.
"""
import pytest

from models.enums import JobStatus
from services.kubejob_service import KubeEventWatcher


class _Condition:
    def __init__(self, type_, status="True"):
        self.type = type_
        self.status = status


class _Job:
    """Minimal stand-in for a V1Job.

    `spec=None` and `status=None` are both reachable in practice -- the API server omits
    status on a Job it has only just accepted -- so they are constructible here rather
    than assumed away.
    """

    def __init__(self, conditions=None, succeeded=None, failed=None,
                 completions=1, has_spec=True, has_status=True):
        self.spec = type("spec", (), {"completions": completions})() if has_spec else None
        if has_status:
            self.status = type("status", (), {
                "conditions": conditions,
                "succeeded": succeeded,
                "failed": failed,
            })()
        else:
            self.status = None


def _derive(**kwargs):
    return KubeEventWatcher._derive_phase(_Job(**kwargs))


class TestSuccessIsDetectedFromTheWholeConditionList:
    def test_complete_alone_is_a_completion(self):
        """NOT DETECTED BEFORE. 'Complete' is the condition every successful Job gets,
        and on clusters older than 1.31 it is the only one. The old test named
        'SuccessCriteriaMet' exclusively, so ordinary successes were missed outright."""
        assert _derive(conditions=[_Condition("Complete")]) == JobStatus.COMPLETED

    def test_success_criteria_met_alone_is_a_completion(self):
        """The one shape the old logic did handle. Kept so widening the check cannot
        quietly drop it."""
        assert _derive(conditions=[_Condition("SuccessCriteriaMet")]) == JobStatus.COMPLETED

    def test_success_criteria_met_then_complete(self):
        """The real 1.31+ ordering: the controller appends SuccessCriteriaMet, then
        Complete."""
        assert _derive(conditions=[
            _Condition("SuccessCriteriaMet"), _Condition("Complete"),
        ]) == JobStatus.COMPLETED

    def test_complete_not_at_index_zero(self):
        """NOT DETECTED BEFORE. A Job that was suspended and later completed carries
        Suspended(False) ahead of Complete. conditions[0] is 'Suspended', so the old
        check missed it and fell through to the failed-counter branch, which -- with
        failed unset -- produced no phase at all."""
        assert _derive(conditions=[
            _Condition("Suspended", status="False"), _Condition("Complete"),
        ]) == JobStatus.COMPLETED


class TestFailureIsDetectedFromTheWholeConditionList:
    def test_failed_condition(self):
        assert _derive(conditions=[_Condition("Failed")]) == JobStatus.ERROR

    def test_failure_target_condition(self):
        """'FailureTarget' is set before 'Failed' while the controller is still tidying
        up pods. Treated as terminal: the outcome is already decided."""
        assert _derive(conditions=[_Condition("FailureTarget")]) == JobStatus.ERROR

    def test_failure_target_then_failed(self):
        assert _derive(conditions=[
            _Condition("FailureTarget"), _Condition("Failed"),
        ]) == JobStatus.ERROR

    def test_deadline_exceeded_with_no_failed_pods(self):
        """NOT DETECTED BEFORE, and the most damaging of the three. activeDeadlineSeconds
        is three days; when it elapses the Job is marked Failed but its pods are deleted
        rather than counted, so `failed` stays 0. The old logic matched neither branch and
        logged '>> Skipped job update', leaving a job that had definitively run out of
        time still reporting 'queued' to the client polling it."""
        assert _derive(
            conditions=[_Condition("Failed")], failed=0,
        ) == JobStatus.ERROR


class TestConditionStatusIsHonored:
    """A condition present with status 'False' asserts the opposite of itself. The old
    code read only `.type` and would have accepted either."""

    def test_complete_false_is_not_a_completion(self):
        assert _derive(conditions=[_Condition("Complete", status="False")]) == JobStatus.PROCESSING

    def test_failed_false_is_not_a_failure(self):
        assert _derive(conditions=[_Condition("Failed", status="False")]) == JobStatus.PROCESSING


class TestNonTerminalShapes:
    def test_no_conditions_yet(self):
        """A freshly accepted Job. Not an error, and not a completion."""
        assert _derive(conditions=None) == JobStatus.PROCESSING

    def test_an_empty_condition_list(self):
        assert _derive(conditions=[]) == JobStatus.PROCESSING

    def test_absent_status(self):
        assert _derive(has_status=False) == JobStatus.PROCESSING

    def test_suspended_is_not_terminal(self):
        """A suspended Job has not finished; it has been paused and may yet run."""
        assert _derive(conditions=[_Condition("Suspended")]) == JobStatus.PROCESSING

    def test_zero_failed_pods_is_not_a_failure(self):
        assert _derive(failed=0) == JobStatus.PROCESSING


class TestCounterFallback:
    """Last resort, for the case where conditions are unavailable or lag behind. Reached
    only when no terminal condition is set."""

    def test_succeeded_meets_a_single_completion(self):
        assert _derive(succeeded=1) == JobStatus.COMPLETED

    def test_succeeded_short_of_a_parallel_completion_target(self):
        """A Job needing 3 successes with 1 recorded is still running. Comparing against
        `completions` rather than testing `succeeded > 0` is what keeps this correct."""
        assert _derive(succeeded=1, completions=3) == JobStatus.PROCESSING

    def test_succeeded_meets_a_parallel_completion_target(self):
        assert _derive(succeeded=3, completions=3) == JobStatus.COMPLETED

    def test_completions_defaults_to_one_when_unset(self):
        """`spec.completions` is None for an ordinary single-pod Job, which must not make
        the comparison fail or throw."""
        assert _derive(succeeded=1, completions=None) == JobStatus.COMPLETED

    def test_completions_defaults_to_one_when_spec_is_absent(self):
        assert _derive(succeeded=1, has_spec=False) == JobStatus.COMPLETED

    def test_a_failed_pod_is_terminal(self):
        """Jobs here are created with backoffLimit=0, so there is no retry to wait for."""
        assert _derive(failed=1) == JobStatus.ERROR


class TestPrecedence:
    def test_success_wins_over_failure(self):
        """Both conditions true is contradictory but observable during a race between the
        controller's writes. Resolved toward completion deliberately: the results of a Job
        that reported success exist and are worth serving, whereas marking it 'error'
        strands output the user could have had."""
        assert _derive(conditions=[
            _Condition("Failed"), _Condition("Complete"),
        ]) == JobStatus.COMPLETED

    def test_a_terminal_condition_wins_over_the_counters(self):
        """DeadlineExceeded again, from the other side: one pod did succeed before the
        deadline killed the Job, but the Job as a whole failed. The condition is
        authoritative; the counters are only a fallback for when it is missing."""
        assert _derive(
            conditions=[_Condition("Failed")], succeeded=1,
        ) == JobStatus.ERROR
