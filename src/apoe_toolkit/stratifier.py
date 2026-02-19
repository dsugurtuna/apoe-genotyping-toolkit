#!/usr/bin/env python3
"""
Cohort Stratifier for APOE-Based Recall Studies
================================================

Generates stratified recall lists for studies that require balanced
representation of APOE genotypes across demographic groups. Supports
complex stratification requirements such as:

  - 50/50 APOE e4 carrier vs. e3/e3 non-carrier split
  - Gender-balanced groups (e.g. female menopausal staging)
  - Age-band matching between case and control arms
  - Exclusion of e2 carriers

This module was inspired by work supporting the NBR267 (Memory and Menopause)
study, which required 816 participants stratified by menopause stage, gender,
and APOE4 carrier status.

Author: Ugur Tuna
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class StratificationConfig:
    """Configuration for a stratified recall list."""

    study_name: str = "Unnamed Study"
    target_female_count: int = 640
    target_male_count: int = 176
    apoe_carrier_ratio: float = 0.5  # proportion of e4 carriers in each group
    exclude_e2_carriers: bool = True
    female_age_bands: Optional[list[tuple[int, int]]] = None
    male_age_bands: Optional[list[tuple[int, int]]] = None

    def __post_init__(self):
        if self.female_age_bands is None:
            self.female_age_bands = [
                (35, 39),
                (40, 44),
                (45, 49),
                (50, 54),
                (55, 59),
                (60, 64),
                (65, 69),
            ]
        if self.male_age_bands is None:
            self.male_age_bands = [
                (35, 39),
                (40, 44),
                (45, 49),
                (50, 54),
                (55, 59),
                (60, 64),
                (65, 69),
            ]


@dataclass
class RecallList:
    """Generated recall list for one arm of the study."""

    arm_name: str
    participants: pd.DataFrame
    summary: dict[str, int] = field(default_factory=dict)


@dataclass
class StratificationResult:
    """Complete output of a stratification run."""

    config: StratificationConfig
    female_list: Optional[RecallList] = None
    male_list: Optional[RecallList] = None
    excluded_count: int = 0
    total_eligible: int = 0

    @property
    def total_selected(self) -> int:
        n = 0
        if self.female_list is not None:
            n += len(self.female_list.participants)
        if self.male_list is not None:
            n += len(self.male_list.participants)
        return n


class CohortStratifier:
    """
    Generates stratified recall lists balanced by APOE status, gender,
    and age band.

    Example::

        stratifier = CohortStratifier()
        config = StratificationConfig(
            study_name="NBR267 Memory and Menopause",
            target_female_count=640,
            target_male_count=176,
            exclude_e2_carriers=True,
        )
        result = stratifier.stratify(cohort_df, config)
    """

    def stratify(
        self,
        cohort: pd.DataFrame,
        config: StratificationConfig,
    ) -> StratificationResult:
        """
        Produce recall lists from a cohort DataFrame.

        Parameters
        ----------
        cohort : pd.DataFrame
            Must contain columns: sample_id, apoe_genotype, gender, age
            (or year_of_birth).
        config : StratificationConfig
            Study-specific stratification parameters.

        Returns
        -------
        StratificationResult
        """
        df = cohort.copy()

        # Normalise column names
        col_map = {c.lower(): c for c in df.columns}
        if "gender" not in col_map and "sex" in col_map:
            df = df.rename(columns={col_map["sex"]: "gender"})

        # Derive age if year_of_birth is present but age is not
        if "age" not in df.columns and "year_of_birth" in df.columns:
            import datetime
            current_year = datetime.datetime.now().year
            df["age"] = current_year - df["year_of_birth"]

        required = {"sample_id", "apoe_genotype", "gender", "age"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Apply exclusions
        initial_count = len(df)
        if config.exclude_e2_carriers:
            df = df[~df["apoe_genotype"].str.contains("e2", case=False, na=False)]

        df = df[df["apoe_genotype"] != "Indeterminate"]
        excluded = initial_count - len(df)

        # Classify carrier status
        df["is_e4_carrier"] = df["apoe_genotype"].str.contains("e4", case=False, na=False)

        result = StratificationResult(
            config=config,
            excluded_count=excluded,
            total_eligible=len(df),
        )

        # Stratify females
        females = df[df["gender"].astype(str).str.upper().isin(["F", "FEMALE", "2"])]
        if not females.empty and config.target_female_count > 0:
            result.female_list = self._select_arm(
                females,
                arm_name="Female",
                target_count=config.target_female_count,
                carrier_ratio=config.apoe_carrier_ratio,
                age_bands=config.female_age_bands,
            )

        # Stratify males
        males = df[df["gender"].astype(str).str.upper().isin(["M", "MALE", "1"])]
        if not males.empty and config.target_male_count > 0:
            result.male_list = self._select_arm(
                males,
                arm_name="Male",
                target_count=config.target_male_count,
                carrier_ratio=config.apoe_carrier_ratio,
                age_bands=config.male_age_bands,
            )

        return result

    def _select_arm(
        self,
        df: pd.DataFrame,
        arm_name: str,
        target_count: int,
        carrier_ratio: float,
        age_bands: Optional[list[tuple[int, int]]],
    ) -> RecallList:
        """Select participants for one arm, balanced by APOE and age."""
        carriers = df[df["is_e4_carrier"]]
        non_carriers = df[~df["is_e4_carrier"]]

        n_carriers = int(target_count * carrier_ratio)
        n_non_carriers = target_count - n_carriers

        if age_bands:
            selected_carriers = self._sample_by_age_bands(carriers, n_carriers, age_bands)
            selected_non = self._sample_by_age_bands(non_carriers, n_non_carriers, age_bands)
        else:
            selected_carriers = carriers.head(min(n_carriers, len(carriers)))
            selected_non = non_carriers.head(min(n_non_carriers, len(non_carriers)))

        combined = pd.concat([selected_carriers, selected_non], ignore_index=True)

        summary = {
            "target": target_count,
            "selected": len(combined),
            "e4_carriers": len(selected_carriers),
            "non_carriers": len(selected_non),
            "available_carriers": len(carriers),
            "available_non_carriers": len(non_carriers),
        }

        return RecallList(arm_name=arm_name, participants=combined, summary=summary)

    @staticmethod
    def _sample_by_age_bands(
        df: pd.DataFrame,
        target: int,
        age_bands: list[tuple[int, int]],
    ) -> pd.DataFrame:
        """Distribute samples proportionally across age bands."""
        per_band = max(1, target // len(age_bands))
        remainder = target - per_band * len(age_bands)

        parts: list[pd.DataFrame] = []
        for i, (lo, hi) in enumerate(age_bands):
            band_df = df[(df["age"] >= lo) & (df["age"] <= hi)]
            n = per_band + (1 if i < remainder else 0)
            parts.append(band_df.head(min(n, len(band_df))))

        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    @staticmethod
    def export_recall_lists(
        result: StratificationResult,
        output_dir: str,
    ) -> list[str]:
        """
        Write recall lists and summary to CSV files.

        Parameters
        ----------
        result : StratificationResult
        output_dir : str

        Returns
        -------
        list of str
            Paths of created files.
        """
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        created: list[str] = []

        study_slug = result.config.study_name.replace(" ", "_").lower()

        if result.female_list is not None:
            p = out / f"{study_slug}_female_recall_list.csv"
            result.female_list.participants.to_csv(p, index=False)
            created.append(str(p))

        if result.male_list is not None:
            p = out / f"{study_slug}_male_recall_list.csv"
            result.male_list.participants.to_csv(p, index=False)
            created.append(str(p))

        # Write summary
        summary_lines = [
            f"Study: {result.config.study_name}",
            f"Total eligible after exclusions: {result.total_eligible:,}",
            f"Excluded (e2/indeterminate): {result.excluded_count:,}",
            f"Total selected: {result.total_selected:,}",
            "",
        ]
        for arm in [result.female_list, result.male_list]:
            if arm is not None:
                summary_lines.append(f"{arm.arm_name} Arm:")
                for k, v in arm.summary.items():
                    summary_lines.append(f"  {k}: {v:,}")
                summary_lines.append("")

        sp = out / f"{study_slug}_recall_summary.txt"
        sp.write_text("\n".join(summary_lines))
        created.append(str(sp))

        return created
