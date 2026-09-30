"""Tests for genotype input handling in APOECaller."""

from pathlib import Path

import pytest

from apoe_toolkit.caller import APOECaller, normalise_genotype

EXAMPLE = Path(__file__).resolve().parent.parent / "data" / "example"


def test_normalise_genotype() -> None:
    assert normalise_genotype("C/T") == "CT"
    assert normalise_genotype("t c") == "CT"
    assert normalise_genotype("TT") == "TT"
    assert normalise_genotype("0 0") == "00"


def test_csv_with_counts_and_sample_id_column() -> None:
    results = APOECaller().call_from_csv(str(EXAMPLE / "example_dosages.csv"))
    calls = {r.sample_id: r.apoe_genotype for r in results}
    # C copies at rs429358, T copies at rs7412.
    assert calls["SYNTH_001"] == "e3/e3"  # 0, 0
    assert calls["SYNTH_002"] == "e2/e2"  # 0, 2
    assert calls["SYNTH_003"] == "e4/e4"  # 2, 0
    assert calls["SYNTH_004"] == "e3/e4"  # 1, 0


def test_csv_with_slash_genotypes(tmp_path: Path) -> None:
    csv = tmp_path / "g.csv"
    csv.write_text("IID,rs429358,rs7412\nS1,T/C,C/C\nS2,T/T,C/T\n")
    calls = [r.apoe_genotype for r in APOECaller().call_from_csv(str(csv))]
    assert calls == ["e3/e4", "e2/e3"]


def test_raw_uses_counted_allele_from_column_name(tmp_path: Path) -> None:
    raw = tmp_path / "x.raw"
    # rs429358 counts T here, not C: 2 copies of T means 0 copies of C.
    raw.write_text(
        "FID IID PAT MAT SEX PHENOTYPE rs429358_T rs7412_T\n"
        "F1 S1 0 0 1 -9 2 0\n"
        "F2 S2 0 0 1 -9 1 0\n"
        "F3 S3 0 0 1 -9 NA 0\n"
    )
    calls = {r.sample_id: r.apoe_genotype for r in APOECaller().call_from_raw(str(raw))}
    assert calls == {"S1": "e3/e3", "S2": "e3/e4", "S3": "Indeterminate"}


def test_raw_on_reverse_strand_is_refused(tmp_path: Path) -> None:
    raw = tmp_path / "x.raw"
    raw.write_text(
        "FID IID PAT MAT SEX PHENOTYPE rs429358_G rs7412_A\nF1 S1 0 0 1 -9 0 0\n"
    )
    with pytest.raises(ValueError, match="reverse strand"):
        APOECaller().call_from_raw(str(raw))
