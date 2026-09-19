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
import io
import pathlib
import re

import yaml

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


def _load_yaml(relative_path):
    """Load a repo YAML file, tolerating the stray trailing tab in chart/values.yaml."""
    text = io.open(REPO_ROOT / relative_path, encoding="utf-8").read()
    return yaml.safe_load("\n".join(line.rstrip("\t") for line in text.split("\n")))


def _simplefold_config(relative_path):
    doc = _load_yaml(relative_path)
    return (doc.get("config") or doc)["kubernetes_jobs"]["ml-simplefold"]


class TestTheVariantIsDeclaredInConfig:
    def test_the_chart_declares_a_model(self):
        """If the key is absent the code falls back to a constant, which is how the
        variant became invisible in the first place. It must be stated where the image
        is stated."""
        assert _simplefold_config("chart/values.yaml").get("simplefoldModel")

    def test_the_in_repo_default_config_declares_a_model(self):
        assert _simplefold_config("app/cfg/config.yaml").get("simplefoldModel")

    def test_both_config_files_agree(self):
        """The chart replaces app/cfg/config.yaml entirely when deployed, so a
        disagreement means local runs and deployed runs use different models -- the same
        class of split this change exists to close."""
        assert (_simplefold_config("chart/values.yaml")["simplefoldModel"]
                == _simplefold_config("app/cfg/config.yaml")["simplefoldModel"])


class TestTheFallbackCannotDriftFromTheDeployedValue:
    def test_the_constant_matches_the_chart(self):
        """`DEFAULT_SIMPLEFOLD_MODEL` is what runs if the config key goes missing. If it
        disagrees with the chart, deleting a line of YAML silently swaps the model --
        exactly the failure this file is about, reintroduced through the back door."""
        source = io.open(REPO_ROOT / "app/services/job_builder.py", encoding="utf-8").read()
        match = re.search(r'DEFAULT_SIMPLEFOLD_MODEL\s*=\s*"([^"]+)"', source)
        assert match, "DEFAULT_SIMPLEFOLD_MODEL not found"
        assert match.group(1) == _simplefold_config("chart/values.yaml")["simplefoldModel"]


class TestTheModelIsNoLongerHardcoded:
    def test_no_literal_variant_remains_in_the_command(self):
        """A second hardcoded variant would re-create the two-sources-of-truth split.
        The command must interpolate, never name a variant directly."""
        source = io.open(REPO_ROOT / "app/services/job_builder.py", encoding="utf-8").read()
        # A real variant name, not an f-string placeholder: `{simplefold_model}` is the
        # fix, `simplefold_100M` would be the bug.
        hardcoded = re.findall(r"--simplefold_model\s+(simplefold_\S+)", source)
        assert hardcoded == [], f"variant hardcoded in command: {hardcoded}"


class TestTheImageTagIsKnownToBeMisleading:
    def test_the_tag_and_the_model_disagree_and_that_is_recorded(self):
        """Not a failure to fix here: the image lives in a personal Docker Hub namespace
        with no public source, so it cannot be retagged from this repo. What matters is
        that nobody reads the tag as the model again, so the config must carry a comment
        saying which one wins. This asserts that explanation is present."""
        cfg = _simplefold_config("chart/values.yaml")
        tag = cfg["image"].rsplit(":", 1)[-1]
        model = cfg["simplefoldModel"]
        tag_digits = re.sub(r"[^0-9.]", "", tag)
        model_digits = re.sub(r"[^0-9.]", "", model.replace("simplefold_", ""))
        if tag_digits and tag_digits != model_digits:
            text = io.open(REPO_ROOT / "chart/values.yaml", encoding="utf-8").read()
            block = text[text.index("ml-simplefold:"):]
            block = block[:block.index("simplefoldModel:")]
            assert "ENTRYPOINT" in block or "entrypoint" in block, (
                "image tag disagrees with the model and nothing explains which wins")
