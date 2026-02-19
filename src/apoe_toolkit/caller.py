#!/usr/bin/env python3
"""
APOE Genotype Caller Module
============================

Implements the standard clinical mapping logic for determining APOE genotypes
from the two defining SNPs: rs429358 and rs7412. Supports multiple input
formats including PLINK .raw, .ped, and pre-extracted CSV.

The APOE gene encodes three major isoforms (E2, E3, E4) defined by two
single-nucleotide polymorphisms on chromosome 19:
  - rs429358 (codon 112): T>C substitution
  - rs7412 (codon 158): C>T substitution

Haplotype definitions:
  - E2: T at rs429358, T at rs7412
  - E3: T at rs429358, C at rs7412
  - E4: C at rs429358, C at rs7412

Author: Ugur Tuna
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Data Models
# ---------------------------------------------------------------------------

APOE_DIPLOTYPE_MAP: dict[tuple[str, str], str] = {
    # Homozygous
    ("TT", "TT"): "e2/e2",
    ("TT", "CC"): "e3/e3",
    ("CC", "CC"): "e4/e4",
    # Heterozygous (common)
    ("TT", "CT"): "e2/e3",
    ("TT", "TC"): "e2/e3",
    ("CT", "CC"): "e3/e4",
    ("TC", "CC"): "e3/e4",
    # Heterozygous (rare / ambiguous — e2/e4 is far more frequent than e1/e3)
    ("CT", "CT"): "e2/e4",
    ("TC", "TC"): "e2/e4",
    ("CT", "TC"): "e2/e4",
    ("TC", "CT"): "e2/e4",
}

RISK_PROFILES: dict[str, str] = {
    "e2/e2": "Reduced risk",
    "e2/e3": "Reduced risk",
    "e3/e3": "Population baseline",
    "e2/e4": "Uncertain / mixed",
    "e3/e4": "Increased risk",
    "e4/e4": "Substantially increased risk",
}


@dataclass
class APOEResult:
    """Represents the APOE genotype call for a single sample."""

    sample_id: str
    rs429358_genotype: str
    rs7412_genotype: str
    apoe_genotype: str
    risk_profile: str
    is_e4_carrier: bool
    is_e2_carrier: bool


@dataclass
class APOESummary:
    """Aggregate statistics for a cohort of APOE calls."""

    total_samples: int = 0
    genotype_counts: dict[str, int] = field(default_factory=dict)
    e4_carrier_count: int = 0
    e2_carrier_count: int = 0
    e3e3_count: int = 0
    indeterminate_count: int = 0


# ---------------------------------------------------------------------------
# Caller Class
# ---------------------------------------------------------------------------

class APOECaller:
    """
    Determines APOE genotypes from genomic data.

    Supports reading from PLINK .raw format (--recode A) as well as
    pre-extracted CSV/TSV files containing the two key SNP columns.

    Example usage::

        caller = APOECaller()
        results = caller.call_from_raw("/path/to/batch28.raw")
        summary = caller.summarise(results)
    """

    REQUIRED_SNPS = ("rs429358", "rs7412")

    def call_from_raw(self, filepath: str) -> list[APOEResult]:
        """
        Call APOE genotypes from a PLINK .raw file.

        The .raw file is expected to be generated via::

            plink --bfile <prefix> --extract snp_list.txt --recode A --out <prefix>

        Parameters
        ----------
        filepath : str
            Path to the PLINK .raw file.

        Returns
        -------
        list of APOEResult
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")

        df = pd.read_csv(path, sep=r"\s+", engine="python")
        return self._call_from_dataframe(df)

    def call_from_ped(self, filepath: str) -> list[APOEResult]:
        """
        Call APOE genotypes from a PLINK .ped file with two extracted SNPs.

        Parameters
        ----------
        filepath : str
            Path to the .ped file (compound genotypes format).

        Returns
        -------
        list of APOEResult
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")

        header = ["FID", "IID", "PAT", "MAT", "SEX", "PHENO", "rs429358", "rs7412"]
        df = pd.read_csv(path, sep=r"\s+", names=header, header=None, engine="python")

        results: list[APOEResult] = []
        for _, row in df.iterrows():
            gt_429 = str(row["rs429358"])
            gt_741 = str(row["rs7412"])
            genotype = APOE_DIPLOTYPE_MAP.get((gt_429, gt_741), "Indeterminate")
            risk = RISK_PROFILES.get(genotype, "Unknown")

            results.append(
                APOEResult(
                    sample_id=str(row["IID"]),
                    rs429358_genotype=gt_429,
                    rs7412_genotype=gt_741,
                    apoe_genotype=genotype,
                    risk_profile=risk,
                    is_e4_carrier="e4" in genotype,
                    is_e2_carrier="e2" in genotype,
                )
            )
        return results

    def call_from_csv(
        self,
        filepath: str,
        sample_col: str = "IID",
        rs429358_col: str = "rs429358",
        rs7412_col: str = "rs7412",
        sep: str = ",",
    ) -> list[APOEResult]:
        """
        Call APOE genotypes from a generic CSV/TSV file.

        Parameters
        ----------
        filepath : str
            Path to the input file.
        sample_col : str
            Column name for sample identifiers.
        rs429358_col : str
            Column name for rs429358 genotype values.
        rs7412_col : str
            Column name for rs7412 genotype values.
        sep : str
            Column delimiter.

        Returns
        -------
        list of APOEResult
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")

        df = pd.read_csv(path, sep=sep)
        required = {sample_col, rs429358_col, rs7412_col}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing columns in input: {missing}")

        results: list[APOEResult] = []
        for _, row in df.iterrows():
            gt_429 = str(row[rs429358_col])
            gt_741 = str(row[rs7412_col])
            genotype = APOE_DIPLOTYPE_MAP.get((gt_429, gt_741), "Indeterminate")
            risk = RISK_PROFILES.get(genotype, "Unknown")

            results.append(
                APOEResult(
                    sample_id=str(row[sample_col]),
                    rs429358_genotype=gt_429,
                    rs7412_genotype=gt_741,
                    apoe_genotype=genotype,
                    risk_profile=risk,
                    is_e4_carrier="e4" in genotype,
                    is_e2_carrier="e2" in genotype,
                )
            )
        return results

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------

    def _call_from_dataframe(self, df: pd.DataFrame) -> list[APOEResult]:
        """Resolve genotype from a DataFrame containing allele-count columns."""
        results: list[APOEResult] = []
        # In .raw format the columns are named <SNP>_<counted_allele>
        snp_cols = [c for c in df.columns if c.startswith(("rs429358", "rs7412"))]
        if len(snp_cols) < 2:
            raise ValueError(
                "Could not find rs429358 and rs7412 columns in the .raw file."
            )

        rs429_col = [c for c in snp_cols if c.startswith("rs429358")][0]
        rs741_col = [c for c in snp_cols if c.startswith("rs7412")][0]

        for _, row in df.iterrows():
            dose_429 = row[rs429_col]
            dose_741 = row[rs741_col]
            genotype, gt_429, gt_741 = self._dosage_to_diplotype(dose_429, dose_741)
            risk = RISK_PROFILES.get(genotype, "Unknown")

            results.append(
                APOEResult(
                    sample_id=str(row.get("IID", row.get("FID", ""))),
                    rs429358_genotype=gt_429,
                    rs7412_genotype=gt_741,
                    apoe_genotype=genotype,
                    risk_profile=risk,
                    is_e4_carrier="e4" in genotype,
                    is_e2_carrier="e2" in genotype,
                )
            )
        return results

    @staticmethod
    def _dosage_to_diplotype(
        dose_429: float, dose_741: float
    ) -> tuple[str, str, str]:
        """
        Convert allele dosages (0/1/2) to diplotype strings and resolve the
        APOE genotype.

        The counted allele in PLINK .raw for rs429358 is C (risk) and for
        rs7412 is T (protective). We convert back to the two-character
        genotype representation used by the lookup table.
        """
        try:
            d429 = int(round(float(dose_429)))
            d741 = int(round(float(dose_741)))
        except (ValueError, TypeError):
            return ("Indeterminate", str(dose_429), str(dose_741))

        gt_429_map = {0: "TT", 1: "CT", 2: "CC"}
        gt_741_map = {0: "CC", 1: "CT", 2: "TT"}

        gt_429 = gt_429_map.get(d429, "??")
        gt_741 = gt_741_map.get(d741, "??")

        genotype = APOE_DIPLOTYPE_MAP.get((gt_429, gt_741), "Indeterminate")
        return genotype, gt_429, gt_741

    # ------------------------------------------------------------------
    # summarisation
    # ------------------------------------------------------------------

    @staticmethod
    def summarise(results: list[APOEResult]) -> APOESummary:
        """
        Produce aggregate statistics from a list of APOE results.

        Parameters
        ----------
        results : list of APOEResult

        Returns
        -------
        APOESummary
        """
        summary = APOESummary(total_samples=len(results))
        for r in results:
            summary.genotype_counts[r.apoe_genotype] = (
                summary.genotype_counts.get(r.apoe_genotype, 0) + 1
            )
            if r.is_e4_carrier:
                summary.e4_carrier_count += 1
            if r.is_e2_carrier:
                summary.e2_carrier_count += 1
            if r.apoe_genotype == "e3/e3":
                summary.e3e3_count += 1
            if r.apoe_genotype == "Indeterminate":
                summary.indeterminate_count += 1
        return summary

    @staticmethod
    def results_to_dataframe(results: list[APOEResult]) -> pd.DataFrame:
        """Convert a list of APOEResult to a pandas DataFrame."""
        records = [
            {
                "sample_id": r.sample_id,
                "rs429358": r.rs429358_genotype,
                "rs7412": r.rs7412_genotype,
                "apoe_genotype": r.apoe_genotype,
                "risk_profile": r.risk_profile,
                "is_e4_carrier": r.is_e4_carrier,
                "is_e2_carrier": r.is_e2_carrier,
            }
            for r in results
        ]
        return pd.DataFrame(records)
