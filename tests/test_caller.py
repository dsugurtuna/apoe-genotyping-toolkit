"""Tests for the APOECaller class."""

import pandas as pd
import pytest

from apoe_toolkit.caller import APOECaller, APOEResult


@pytest.fixture
def caller():
    return APOECaller()


class TestDosageToDiplotype:
    """Test the internal dosage-to-diplotype conversion (returns a tuple)."""

    def test_e3_e3(self, caller):
        """rs429358=0 (TT), rs7412=0 (CC) => e3/e3."""
        genotype, gt_429, gt_741 = caller._dosage_to_diplotype(0, 0)
        assert genotype == "e3/e3"

    def test_e2_e2(self, caller):
        genotype, _, _ = caller._dosage_to_diplotype(0, 2)
        assert genotype == "e2/e2"

    def test_e4_e4(self, caller):
        genotype, _, _ = caller._dosage_to_diplotype(2, 0)
        assert genotype == "e4/e4"

    def test_e3_e4(self, caller):
        genotype, _, _ = caller._dosage_to_diplotype(1, 0)
        assert genotype == "e3/e4"

    def test_e2_e3(self, caller):
        genotype, _, _ = caller._dosage_to_diplotype(0, 1)
        assert genotype == "e2/e3"

    def test_e2_e4(self, caller):
        genotype, _, _ = caller._dosage_to_diplotype(1, 1)
        assert genotype == "e2/e4"

    def test_indeterminate_on_invalid(self, caller):
        genotype, _, _ = caller._dosage_to_diplotype("bad", "data")
        assert genotype == "Indeterminate"


class TestCallFromCSV:
    """Test the CSV ingestion path using the IID column default."""

    def test_basic_csv(self, caller, tmp_path):
        csv = tmp_path / "test.csv"
        csv.write_text(
            "IID,rs429358,rs7412\n"
            "S1,TT,CC\n"
            "S2,CC,CC\n"
            "S3,TT,TT\n"
        )
        results = caller.call_from_csv(str(csv))
        assert len(results) == 3
        assert results[0].apoe_genotype == "e3/e3"
        assert results[1].apoe_genotype == "e4/e4"
        assert results[2].apoe_genotype == "e2/e2"

    def test_custom_columns(self, caller, tmp_path):
        csv = tmp_path / "test.csv"
        csv.write_text(
            "sample_id,rs429358,rs7412\n"
            "S1,TT,CC\n"
        )
        results = caller.call_from_csv(str(csv), sample_col="sample_id")
        assert len(results) == 1
        assert results[0].sample_id == "S1"


class TestSummarise:
    def test_summary_counts(self, caller, tmp_path):
        csv = tmp_path / "test.csv"
        csv.write_text(
            "IID,rs429358,rs7412\n"
            "S1,TT,CC\n"
            "S2,TT,CC\n"
            "S3,CC,CC\n"
        )
        results = caller.call_from_csv(str(csv))
        summary = caller.summarise(results)
        assert summary.total_samples == 3
        assert summary.genotype_counts["e3/e3"] == 2
        assert summary.genotype_counts["e4/e4"] == 1
        assert summary.e4_carrier_count == 1


class TestResultsToDataFrame:
    def test_dataframe_output(self, caller, tmp_path):
        csv = tmp_path / "test.csv"
        csv.write_text("IID,rs429358,rs7412\nS1,TT,CC\n")
        results = caller.call_from_csv(str(csv))
        df = caller.results_to_dataframe(results)
        assert isinstance(df, pd.DataFrame)
        assert "sample_id" in df.columns
        assert "apoe_genotype" in df.columns
        assert len(df) == 1
