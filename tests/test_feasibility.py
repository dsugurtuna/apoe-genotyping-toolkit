"""Tests for the APOEFeasibilityEstimator class."""

import pytest

from apoe_toolkit.caller import APOECaller, APOEResult
from apoe_toolkit.feasibility import APOEFeasibilityEstimator


@pytest.fixture
def estimator():
    return APOEFeasibilityEstimator()


@pytest.fixture
def sample_results():
    """Build a small set of APOE results for testing."""
    data = [
        ("S1", "TT", "CC", "e3/e3", "Population baseline", False, False),
        ("S2", "CC", "CC", "e4/e4", "Substantially increased risk", True, False),
        ("S3", "CT", "CC", "e3/e4", "Increased risk", True, False),
        ("S4", "TT", "CT", "e2/e3", "Reduced risk", False, True),
        ("S5", "CT", "CT", "e2/e4", "Uncertain / mixed", True, True),
        ("S6", "TT", "CC", "e3/e3", "Population baseline", False, False),
        ("S7", "CT", "CC", "e3/e4", "Increased risk", True, False),
        ("S8", "TT", "CC", "e3/e3", "Population baseline", False, False),
        ("S9", "CC", "CC", "e4/e4", "Substantially increased risk", True, False),
        ("S10", "TT", "CC", "e3/e3", "Population baseline", False, False),
    ]
    return [
        APOEResult(
            sample_id=d[0],
            rs429358_genotype=d[1],
            rs7412_genotype=d[2],
            apoe_genotype=d[3],
            risk_profile=d[4],
            is_e4_carrier=d[5],
            is_e2_carrier=d[6],
        )
        for d in data
    ]


class TestFeasibilityEstimator:
    def test_estimate_from_results(self, estimator, sample_results):
        report = estimator.estimate_from_results(
            results=sample_results,
            target_genotypes=["e4/e4", "e3/e4"],
        )
        assert report.total_genotyped == 10
        assert report.eligible_count == 4  # 2 e4/e4 + 2 e3/e4
        assert 0.0 <= report.eligibility_rate <= 1.0

    def test_estimate_from_csv(self, estimator, tmp_path):
        csv = tmp_path / "results.csv"
        csv.write_text(
            "IID,rs429358,rs7412\n"
            "S1,TT,CC\n"
            "S2,CC,CC\n"
            "S3,CT,CC\n"
            "S4,TT,CT\n"
            "S5,TT,CC\n"
        )
        report = estimator.estimate_from_csv(
            filepath=str(csv),
            target_genotypes=["e4/e4"],
        )
        assert report.total_genotyped == 5
        assert report.eligible_count == 1

    def test_format_report(self, estimator, sample_results):
        report = estimator.estimate_from_results(
            results=sample_results,
            target_genotypes=["e4/e4"],
        )
        text = estimator.format_report(report)
        assert "e4/e4" in text
        assert "Total genotyped" in text

    def test_empty_results(self, estimator):
        report = estimator.estimate_from_results(
            results=[],
            target_genotypes=["e4/e4"],
        )
        assert report.total_genotyped == 0
        assert report.eligible_count == 0
