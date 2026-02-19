#!/usr/bin/env python3
"""
APOE Feasibility Estimator
===========================

Provides rapid feasibility estimates for clinical trials requiring APOE-stratified
participants. Given a genotyped cohort, calculates the number of available
participants meeting specific APOE genotype criteria (e.g. e4/e4 homozygotes
or e3/e4 heterozygotes for an Alzheimer's disease trial).

This module was inspired by work supporting pharmaceutical feasibility enquiries
such as those from NewAmsterdam Pharma (Alzheimer's Disease clinical trial
screening), where approximate counts of e4 carriers were required to assess
study viability.

Author: Ugur Tuna
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

import pandas as pd

from apoe_toolkit.caller import APOECaller, APOEResult

logger = logging.getLogger(__name__)


@dataclass
class FeasibilityReport:
    """Summary of APOE feasibility analysis for a study."""

    study_name: str
    total_genotyped: int = 0
    eligible_count: int = 0
    target_genotypes: List[str] = field(default_factory=list)
    genotype_breakdown: Dict[str, int] = field(default_factory=dict)
    exclusion_criteria: List[str] = field(default_factory=list)
    excluded_count: int = 0
    notes: str = ""

    @property
    def eligibility_rate(self) -> float:
        """Proportion of genotyped cohort that meets criteria."""
        if self.total_genotyped == 0:
            return 0.0
        return self.eligible_count / self.total_genotyped

    def to_dict(self) -> dict:
        """Convert to a plain dictionary for serialisation."""
        return {
            "study_name": self.study_name,
            "total_genotyped": self.total_genotyped,
            "eligible_count": self.eligible_count,
            "eligibility_rate": round(self.eligibility_rate, 4),
            "target_genotypes": self.target_genotypes,
            "genotype_breakdown": self.genotype_breakdown,
            "exclusion_criteria": self.exclusion_criteria,
            "excluded_count": self.excluded_count,
            "notes": self.notes,
        }


class APOEFeasibilityEstimator:
    """
    Estimates participant availability for APOE-stratified clinical studies.

    Typical workflow::

        estimator = APOEFeasibilityEstimator()

        # For an Alzheimer's trial requiring e4 carriers
        report = estimator.estimate_from_results(
            results=apoe_results,
            study_name="NewAmsterdam AD Trial",
            target_genotypes=["e3/e4", "e4/e4"],
            exclude_genotypes=["e2/e2", "e2/e3", "e2/e4"],
        )
        print(report.eligible_count)
    """

    def estimate_from_results(
        self,
        results: List[APOEResult],
        study_name: str = "Unnamed Study",
        target_genotypes: Optional[List[str]] = None,
        exclude_genotypes: Optional[List[str]] = None,
        exclude_indeterminate: bool = True,
    ) -> FeasibilityReport:
        """
        Estimate feasibility from a list of pre-called APOE results.

        Parameters
        ----------
        results : list of APOEResult
            Pre-called APOE genotype results.
        study_name : str
            Label for the feasibility report.
        target_genotypes : list of str, optional
            Genotypes to count as eligible (e.g. ["e3/e4", "e4/e4"]).
            If None, all non-excluded genotypes are counted.
        exclude_genotypes : list of str, optional
            Genotypes to exclude from eligibility.
        exclude_indeterminate : bool
            Whether to exclude samples with indeterminate calls.

        Returns
        -------
        FeasibilityReport
        """
        target_set: Optional[Set[str]] = (
            set(target_genotypes) if target_genotypes else None
        )
        exclude_set: Set[str] = set(exclude_genotypes) if exclude_genotypes else set()
        if exclude_indeterminate:
            exclude_set.add("Indeterminate")

        report = FeasibilityReport(
            study_name=study_name,
            total_genotyped=len(results),
            target_genotypes=target_genotypes or [],
            exclusion_criteria=list(exclude_set),
        )

        eligible_count = 0
        excluded_count = 0

        for r in results:
            gt = r.apoe_genotype
            report.genotype_breakdown[gt] = report.genotype_breakdown.get(gt, 0) + 1

            if gt in exclude_set:
                excluded_count += 1
                continue

            if target_set is None or gt in target_set:
                eligible_count += 1

        report.eligible_count = eligible_count
        report.excluded_count = excluded_count
        return report

    def estimate_from_csv(
        self,
        filepath: str,
        study_name: str = "Unnamed Study",
        target_genotypes: Optional[List[str]] = None,
        exclude_genotypes: Optional[List[str]] = None,
        sample_col: str = "IID",
        rs429358_col: str = "rs429358",
        rs7412_col: str = "rs7412",
        sep: str = ",",
    ) -> FeasibilityReport:
        """
        Run end-to-end feasibility estimation from a CSV file.

        Parameters
        ----------
        filepath : str
            Path to a CSV containing genotype columns.
        study_name : str
            Label for the report.
        target_genotypes : list of str, optional
            Genotypes considered eligible.
        exclude_genotypes : list of str, optional
            Genotypes to exclude.
        sample_col, rs429358_col, rs7412_col, sep
            Column mappings and delimiter.

        Returns
        -------
        FeasibilityReport
        """
        caller = APOECaller()
        results = caller.call_from_csv(
            filepath,
            sample_col=sample_col,
            rs429358_col=rs429358_col,
            rs7412_col=rs7412_col,
            sep=sep,
        )
        return self.estimate_from_results(
            results,
            study_name=study_name,
            target_genotypes=target_genotypes,
            exclude_genotypes=exclude_genotypes,
        )

    @staticmethod
    def format_report(report: FeasibilityReport) -> str:
        """
        Return a human-readable text summary suitable for email or Jira.

        Parameters
        ----------
        report : FeasibilityReport

        Returns
        -------
        str
        """
        lines = [
            f"APOE Feasibility Report: {report.study_name}",
            "=" * 60,
            f"Total genotyped samples  : {report.total_genotyped:,}",
            f"Eligible participants    : {report.eligible_count:,}",
            f"Eligibility rate         : {report.eligibility_rate:.1%}",
            f"Excluded participants    : {report.excluded_count:,}",
            "",
            "Target genotypes         : "
            + (", ".join(report.target_genotypes) if report.target_genotypes else "All"),
            "Exclusion criteria       : "
            + (", ".join(report.exclusion_criteria) if report.exclusion_criteria else "None"),
            "",
            "Genotype breakdown:",
            "-" * 40,
        ]
        for gt, count in sorted(
            report.genotype_breakdown.items(), key=lambda x: -x[1]
        ):
            pct = count / report.total_genotyped * 100 if report.total_genotyped else 0
            lines.append(f"  {gt:15s}  {count:>8,}  ({pct:5.1f}%)")

        if report.notes:
            lines.extend(["", f"Notes: {report.notes}"])

        return "\n".join(lines)
