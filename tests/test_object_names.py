"""Caller-chosen parts of `{job_id}/in/{filename}` (#157 review, round 3)."""
import pytest

from services.object_names import is_plain_file_name, is_valid_job_id


@pytest.mark.parametrize("job_id", [
    "4b3a85ee4efa42c3b5e5c1641e0588f4",   # uuid4 hex: the server, and crispr-copies' randomUUID
    "1793160805123.4417382209",           # crispr-copies' fallback, Date.now() + Math.random()
    "my-job", "sub_i",
])
def test_ids_in_use_are_valid(job_id):
    assert is_valid_job_id(job_id)


@pytest.mark.parametrize("job_id", ["victim/in/x", "victim/in/..", "..", ".", "a\\b", "", "x" * 65, "a b"])
def test_ids_that_could_reach_another_jobs_prefix_are_not(job_id):
    assert not is_valid_job_id(job_id)


@pytest.mark.parametrize("name", ["mutations.csv", "a.b.pdf", "..hidden"])
def test_plain_names(name):
    assert is_plain_file_name(name)


@pytest.mark.parametrize("name", ["../in/f", "x/../job.json", "a\\b", ".", "..", ""])
def test_names_with_a_path_are_not_plain(name):
    assert not is_plain_file_name(name)
