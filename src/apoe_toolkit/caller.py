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
import re
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


def normalise_genotype(value: object) -> str:
    """Return a two-allele genotype as sorted letters: 'C/T', 'T C' -> 'CT'.

    Values that are not two alleles (missing data, '0 0') are returned
    upper-cased without separators and will not match the lookup table.
    """
    text = re.sub(r"[\s/|]", "", str(value)).upper()
    return "".join(sorted(text)) if len(text) == 2 else text


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
            gt_429 = normalise_genotype(row["rs429358"])
            gt_741 = normalise_genotype(row["rs7412"])
            genotype = APOE_DIPLOTYPE_MAP.get((gt_429, gt_741), "Indeterminate")
            results.append(self._result(str(row["IID"]), gt_429, gt_741, genotype))
        return results

    def call_from_csv(
        self,
        filepath: str,
        sample_col: str | None = None,
        rs429358_col: str = "rs429358",
        rs7412_col: str = "rs7412",
        sep: str = ",",
    ) -> list[APOEResult]:
        """
        Call APOE genotypes from a CSV/TSV file.

        The two SNP columns may hold genotype strings (``TT``, ``C/T``, ...)
        or allele counts (0/1/2). Counts are read as the number of C alleles
        at rs429358 and of T alleles at rs7412, the alleles that define e4
        and e2.

        Parameters
        ----------
        filepath : str
            Path to the input file.
        sample_col : str, optional
            Sample ID column. Defaults to ``IID``, or ``sample_id`` if there
            is no ``IID`` column.
        rs429358_col, rs7412_col : str
            Column names for the two SNPs.
        sep : str
            Column delimiter.
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")

        df = pd.read_csv(path, sep=sep)
        if sample_col is None:
            sample_col = "IID" if "IID" in df.columns else "sample_id"
        required = {sample_col, rs429358_col, rs7412_col}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing columns in input: {sorted(missing)}")

        as_counts = pd.api.types.is_numeric_dtype(
            df[rs429358_col]
        ) and pd.api.types.is_numeric_dtype(df[rs7412_col])

        results: list[APOEResult] = []
        for _, row in df.iterrows():
            if as_counts:
                genotype, gt_429, gt_741 = self._dosage_to_diplotype(
                    row[rs429358_col], row[rs7412_col]
                )
            else:
                gt_429 = normalise_genotype(row[rs429358_col])
                gt_741 = normalise_genotype(row[rs7412_col])
                genotype = APOE_DIPLOTYPE_MAP.get((gt_429, gt_741), "Indeterminate")
            results.append(self._result(str(row[sample_col]), gt_429, gt_741, genotype))
        return results

    @staticmethod
    def _result(sample_id: str, gt_429: str, gt_741: str, genotype: str) -> APOEResult:
        return APOEResult(
            sample_id=sample_id,
            rs429358_genotype=gt_429,
            rs7412_genotype=gt_741,
            apoe_genotype=genotype,
            risk_profile=RISK_PROFILES.get(genotype, "Unknown"),
            is_e4_carrier="e4" in genotype,
            is_e2_carrier="e2" in genotype,
        )

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _find_raw_column(columns: list[str], snp: str) -> tuple[str, str | None]:
        """Return the .raw column for ``snp`` and its counted allele, if named."""
        for col in columns:
            if col == snp:
                return col, None
            if col.startswith(snp + "_"):
                return col, col[len(snp) + 1 :].upper()
        raise ValueError(f"Could not find a {snp} column in the .raw file.")

    def _call_from_dataframe(self, df: pd.DataFrame) -> list[APOEResult]:
        """Resolve genotypes from a PLINK .raw DataFrame.

        PLINK names each column ``<SNP>_<counted allele>`` and counts its A1
        allele, which is usually but not always the minor one. The counted
        allele is read from the column name and the count converted to C
        copies at rs429358 and T copies at rs7412. Alleles other than C and T
        usually mean the data are on the reverse strand (A/G), which this
        caller refuses rather than guesses.
        """
        columns = [str(c) for c in df.columns]
        col_429, allele_429 = self._find_raw_column(columns, "rs429358")
        col_741, allele_741 = self._find_raw_column(columns, "rs7412")
        for snp, allele in (("rs429358", allele_429), ("rs7412", allele_741)):
            if allele not in (None, "C", "T"):
                raise ValueError(
                    f"{snp} counts allele {allele!r}; expected C or T. The data "
                    "may be on the reverse strand (A/G) and need flipping first."
                )

        results: list[APOEResult] = []
        for _, row in df.iterrows():
            d429, d741 = row[col_429], row[col_741]
            if allele_429 == "T" and pd.notna(d429):
                d429 = 2 - float(d429)
            if allele_741 == "C" and pd.notna(d741):
                d741 = 2 - float(d741)
            genotype, gt_429, gt_741 = self._dosage_to_diplotype(d429, d741)
            sample_id = str(row["IID"] if "IID" in row else row.get("FID", ""))
            results.append(self._result(sample_id, gt_429, gt_741, genotype))
        return results

    @staticmethod
    def _dosage_to_diplotype(
        dose_429: object, dose_741: object
    ) -> tuple[str, str, str]:
        """
        Convert allele counts to genotype strings and the APOE genotype.

        ``dose_429`` is the number of C alleles at rs429358 and ``dose_741``
        the number of T alleles at rs7412 (0, 1 or 2; non-integers are
        rounded). Missing or out-of-range values give "Indeterminate".
        """
        try:
            d429 = round(float(str(dose_429)))
            d741 = round(float(str(dose_741)))
        except ValueError:
            return ("Indeterminate", str(dose_429), str(dose_741))

        gt_429 = {0: "TT", 1: "CT", 2: "CC"}.get(d429, "??")
        gt_741 = {0: "CC", 1: "CT", 2: "TT"}.get(d741, "??")
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
