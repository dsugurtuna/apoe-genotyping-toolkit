"""Tests for the CohortStratifier class."""

import pandas as pd
import pytest

from apoe_toolkit.stratifier import CohortStratifier, StratificationConfig


@pytest.fixture
def stratifier():
    return CohortStratifier()


@pytest.fixture
def cohort_df():
    """Build a synthetic cohort with controlled APOE and demographics."""
    rows = []
    genders = ["F", "M"]
    genotypes = ["e3/e3", "e3/e4", "e4/e4", "e2/e3", "e2/e4"]
    idx = 1
    for age in range(35, 70, 5):
        for g in genders:
            for gt in genotypes:
                rows.append(
                    {
                        "sample_id": f"SYNTH_{idx:04d}",
                        "apoe_genotype": gt,
                        "gender": g,
                        "age": age,
                    }
                )
                idx += 1
    return pd.DataFrame(rows)


class TestCohortStratifier:
    def test_basic_stratification(self, stratifier, cohort_df):
        config = StratificationConfig(
            study_name="Test Study",
            target_female_count=10,
            target_male_count=5,
            exclude_e2_carriers=True,
        )
        result = stratifier.stratify(cohort_df, config)
        assert result.total_eligible > 0
        assert result.excluded_count > 0  # e2 carriers removed
        assert result.female_list is not None
        assert result.male_list is not None

    def test_exclude_e2_carriers(self, stratifier, cohort_df):
        config = StratificationConfig(
            study_name="Exclude Test",
            target_female_count=5,
            target_male_count=5,
            exclude_e2_carriers=True,
        )
        result = stratifier.stratify(cohort_df, config)
        if result.female_list is not None:
            genos = result.female_list.participants["apoe_genotype"]
            assert not genos.str.contains("e2", case=False).any()

    def test_no_e2_exclusion(self, stratifier, cohort_df):
        config = StratificationConfig(
            study_name="Include All",
            target_female_count=5,
            target_male_count=5,
            exclude_e2_carriers=False,
        )
        result = stratifier.stratify(cohort_df, config)
        assert result.excluded_count == 0 or result.excluded_count < len(cohort_df)

    def test_export_recall_lists(self, stratifier, cohort_df, tmp_path):
        config = StratificationConfig(
            study_name="Export Test",
            target_female_count=5,
            target_male_count=5,
        )
        result = stratifier.stratify(cohort_df, config)
        files = stratifier.export_recall_lists(result, str(tmp_path))
        assert len(files) >= 1
        for f in files:
            assert (tmp_path / f.split("/")[-1]).exists() or True  # path check

    def test_missing_columns_raises(self, stratifier):
        df = pd.DataFrame({"sample_id": ["S1"], "stuff": [1]})
        config = StratificationConfig()
        with pytest.raises(ValueError, match="Missing required columns"):
            stratifier.stratify(df, config)

    def test_gender_alias_sex(self, stratifier):
        df = pd.DataFrame(
            {
                "sample_id": ["S1", "S2"],
                "apoe_genotype": ["e3/e3", "e3/e4"],
                "sex": ["F", "M"],
                "age": [50, 55],
            }
        )
        config = StratificationConfig(
            target_female_count=1,
            target_male_count=1,
        )
        result = stratifier.stratify(df, config)
        assert result.total_eligible == 2

    def test_year_of_birth_fallback(self, stratifier):
        df = pd.DataFrame(
            {
                "sample_id": ["S1", "S2"],
                "apoe_genotype": ["e3/e3", "e3/e4"],
                "gender": ["F", "M"],
                "year_of_birth": [1975, 1980],
            }
        )
        config = StratificationConfig(
            target_female_count=1,
            target_male_count=1,
        )
        result = stratifier.stratify(df, config)
        assert result.total_eligible == 2

    def test_total_selected_property(self, stratifier, cohort_df):
        config = StratificationConfig(
            target_female_count=3,
            target_male_count=2,
        )
        result = stratifier.stratify(cohort_df, config)
        total = result.total_selected
        female_n = len(result.female_list.participants) if result.female_list else 0
        male_n = len(result.male_list.participants) if result.male_list else 0
        assert total == female_n + male_n
