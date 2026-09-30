"""Tests for the feasibility interval, age-band sampling and the CLI."""

from pathlib import Path

import pandas as pd
import pytest

from apoe_toolkit.cli import main
from apoe_toolkit.feasibility import FeasibilityReport
from apoe_toolkit.stratifier import CohortStratifier

EXAMPLE = Path(__file__).resolve().parent.parent / "data" / "example"


def test_wilson_interval() -> None:
    report = FeasibilityReport(study_name="x", total_genotyped=100, eligible_count=50)
    low, high = report.eligibility_interval(0.95)
    assert low == pytest.approx(0.4038, abs=1e-3)
    assert high == pytest.approx(0.5962, abs=1e-3)
    assert FeasibilityReport(study_name="x").eligibility_interval() == (0.0, 0.0)


def test_age_band_sampling_never_exceeds_target() -> None:
    df = pd.DataFrame({"age": [37, 42, 47, 52, 57, 62, 67] * 3})
    bands = [(35, 39), (40, 44), (45, 49), (50, 54), (55, 59), (60, 64), (65, 69)]
    assert len(CohortStratifier._sample_by_age_bands(df, 3, bands)) == 3
    assert len(CohortStratifier._sample_by_age_bands(df, 0, bands)) == 0
    assert len(CohortStratifier._sample_by_age_bands(df, 10, bands)) == 10


def test_cli_feasibility_runs(capsys: pytest.CaptureFixture[str]) -> None:
    main(["feasibility", "--input", str(EXAMPLE / "example_dosages.csv")])
    out = capsys.readouterr().out
    assert "Eligible participants" in out
    assert "95% Wilson interval" in out


def test_cli_call_then_feasibility(tmp_path: Path) -> None:
    out = tmp_path / "calls.csv"
    main(["call", "--input", str(EXAMPLE / "example_dosages.csv"), "-o", str(out)])
    # The output of 'call' (sample_id plus genotype strings) is valid input.
    main(["feasibility", "--input", str(out), "--targets", "e4/e4"])


def test_cli_can_keep_e2_carriers(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    args = ["stratify", "--input", str(EXAMPLE / "example_cohort.csv"),
            "--females", "20", "--males", "10", "-o", str(tmp_path)]  # fmt: skip
    main(args)
    excluded = capsys.readouterr().out
    main([*args, "--no-exclude-e2"])
    kept = capsys.readouterr().out
    assert "Excluded: 0" not in excluded
    assert "Excluded: 0" in kept
