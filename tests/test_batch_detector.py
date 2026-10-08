import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "examples"))

import batch_detector as bd  # noqa: E402
from make_example_data import make_data  # noqa: E402

EXAMPLES = ROOT / "examples"

# small settings so the whole test file runs in a few seconds
FAST = ["--n-perms", "199", "--power-reps", "10", "--power-perms-inner", "19",
        "--calib-reps", "4", "--calib-iters", "8", "--n-grid", "60,120"]


def run_tool(out_dir, *extra):
    bd.main([
        "--input", str(EXAMPLES / "example_features.tsv"),
        "--labels", str(EXAMPLES / "example_metadata.tsv"),
        "--batch", "cohort", "--clr",
        "--output", str(out_dir),
        *FAST, *extra,
    ])
    with open(Path(out_dir) / "batch_detector_report.json") as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    return run_tool(tmp_path_factory.mktemp("default_run"))


def test_chance_expectation():
    assert bd.chance_expectation_r2(2, 39) == pytest.approx(1 / 38)
    assert bd.chance_expectation_r2(4, 283) == pytest.approx(3 / 282)
    assert bd.chance_expectation_r2(3, 120) == pytest.approx(2 / 119)


def test_chance_level_matches_shuffled_labels():
    # the mean R2 over shuffled labels should sit at (g - 1) / (n - 1)
    rng = np.random.default_rng(0)
    n, g = 60, 3
    X = rng.normal(size=(n, 10))
    D = bd.aitchison_dist(X)
    labels = np.repeat(np.arange(g), n // g)
    r2 = [bd.permanova(D, rng.permutation(labels), n_perms=1)["R2"] for _ in range(300)]
    assert np.mean(r2) == pytest.approx(bd.chance_expectation_r2(g, n), abs=0.005)


def test_permanova_separated_groups():
    rng = np.random.default_rng(1)
    X = np.vstack([rng.normal(0, 1, (20, 5)), rng.normal(6, 1, (20, 5))])
    grp = np.array([0] * 20 + [1] * 20)
    res = bd.permanova(bd.aitchison_dist(X), grp, n_perms=199)
    assert res["R2"] > 0.5
    assert res["p_value"] < 0.05


def test_example_data_is_reproducible():
    X, meta = make_data()
    X_file = pd.read_csv(EXAMPLES / "example_features.tsv", sep="\t", index_col=0)
    meta_file = pd.read_csv(EXAMPLES / "example_metadata.tsv", sep="\t", index_col=0)
    pd.testing.assert_frame_equal(X, X_file, check_dtype=False)
    pd.testing.assert_frame_equal(meta, meta_file)


def test_example_cohort_effect(report):
    assert report["n_samples"] == 120
    assert report["n_batch_groups"] == 3
    assert report["batch_delta_R2"] > 0.1
    assert report["permanova_p_batch"] < 0.05


def test_example_no_response_effect(report):
    assert abs(report["response_delta_R2"]) < 0.02
    assert report["permanova_p_response"] > 0.05


def test_example_is_no_go(report):
    assert report["pooling_recommendation"] == "NO-GO"
    assert report["recommendation_version"] == "v2_chance_calibrated"


def test_report_has_expected_fields(report):
    for key in ["response_R2", "response_E0", "response_delta_R2", "batch_R2",
                "batch_delta_R2", "batch_signal_ratio", "joint_margin_model",
                "rho_permutation_null", "pooling_recommendation",
                "recommendation_reason", "min_n_for_80pct_power"]:
        assert key in report


def test_multiplicative_replacement_runs(tmp_path):
    rep = run_tool(tmp_path, "--zero-handling", "multiplicative_replacement",
                   "--pseudocount", "1e-5")
    assert rep["zero_handling"] == "multiplicative_replacement"
    assert rep["batch_delta_R2"] > 0.1


def test_single_cohort_mode(tmp_path):
    bd.main([
        "--input", str(EXAMPLES / "example_features.tsv"),
        "--labels", str(EXAMPLES / "example_metadata.tsv"),
        "--clr", "--output", str(tmp_path), *FAST,
    ])
    with open(tmp_path / "batch_detector_report.json") as fh:
        rep = json.load(fh)
    assert rep["batch_R2"] is None
    assert rep["n_batch_groups"] is None


@pytest.mark.parametrize("rho, pct, expected", [
    (0.1, 40.0, "GO"),
    (1.0, 40.0, "CAUTION"),
    (5.0, 40.0, "NO-GO"),
    (1.0, 99.0, "NO-GO"),
])
def test_recommendation_rule(rho, pct, expected):
    call, reason = bd._pooling_recommendation_v2(
        response_delta_r2=0.02, response_p=0.01, rho=rho, rho_percentile=pct)
    assert call == expected
    assert "simulation-calibrated" not in reason


def test_recommendation_no_excess_signal():
    call, _ = bd._pooling_recommendation_v2(-0.001, 0.4, 3.0, 50.0)
    assert call == "NO-GO"
    call, _ = bd._pooling_recommendation_v2(0.003, 0.2, 3.0, 50.0)
    assert call == "NO-GO"


def test_help_and_version():
    out = subprocess.run([sys.executable, str(ROOT / "batch_detector.py"), "--help"],
                         capture_output=True, text=True)
    assert out.returncode == 0
    assert "--batch" in out.stdout and "--legacy-power" in out.stdout
    out = subprocess.run([sys.executable, str(ROOT / "batch_detector.py"), "--version"],
                         capture_output=True, text=True)
    assert out.stdout.strip().endswith(bd.__version__)
    assert bd.__version__ == "2.0.0"


def test_no_hardcoded_paths():
    text = (ROOT / "batch_detector.py").read_text()
    assert "/Users/" not in text and "/home/" not in text
