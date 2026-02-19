#!/usr/bin/env python3
"""
APOE Genotyping Toolkit — Command-Line Interface
=================================================

Provides three sub-commands:

  call       Determine APOE diplotypes from genotype data.
  feasibility  Estimate the number of participants matching
               target diplotypes for clinical-trial screening.
  stratify   Generate recall lists balanced by APOE status,
             gender, and age band.

Usage::

    apoe-toolkit call  --input data.csv --format csv --output results.csv
    apoe-toolkit feasibility --input results.csv --targets e4/e4 e3/e4
    apoe-toolkit stratify --input cohort.csv --study "NBR267" --females 640 --males 176

Author: Ugur Tuna
"""

import argparse
import sys

from apoe_toolkit.caller import APOECaller
from apoe_toolkit.feasibility import APOEFeasibilityEstimator
from apoe_toolkit.stratifier import CohortStratifier, StratificationConfig


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="apoe-toolkit",
        description="APOE Genotyping Toolkit — call, estimate feasibility, and stratify cohorts.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ---- call ----------------------------------------------------------------
    call_p = sub.add_parser("call", help="Determine APOE diplotypes from input data.")
    call_p.add_argument("--input", "-i", required=True, help="Path to input file.")
    call_p.add_argument(
        "--format",
        "-f",
        choices=["csv", "ped", "raw"],
        default="csv",
        help="Input file format (default: csv).",
    )
    call_p.add_argument("--output", "-o", default=None, help="Path to output CSV.")
    call_p.add_argument(
        "--summary",
        action="store_true",
        help="Print a summary of diplotype counts.",
    )

    # ---- feasibility ---------------------------------------------------------
    feas_p = sub.add_parser(
        "feasibility",
        help="Estimate participant counts for target diplotypes.",
    )
    feas_p.add_argument("--input", "-i", required=True, help="Path to results CSV.")
    feas_p.add_argument(
        "--targets",
        nargs="+",
        default=["e4/e4", "e3/e4"],
        help="Target diplotypes (default: e4/e4 e3/e4).",
    )
    feas_p.add_argument(
        "--confidence",
        type=float,
        default=0.95,
        help="Confidence level for interval (default: 0.95).",
    )

    # ---- stratify ------------------------------------------------------------
    strat_p = sub.add_parser(
        "stratify",
        help="Generate recall lists stratified by APOE, gender, age.",
    )
    strat_p.add_argument("--input", "-i", required=True, help="Path to cohort CSV.")
    strat_p.add_argument("--study", default="Unnamed Study", help="Study name.")
    strat_p.add_argument(
        "--females", type=int, default=640, help="Target female count."
    )
    strat_p.add_argument("--males", type=int, default=176, help="Target male count.")
    strat_p.add_argument(
        "--output-dir",
        "-o",
        default="recall_output",
        help="Output directory for recall lists.",
    )
    strat_p.add_argument(
        "--exclude-e2",
        action="store_true",
        default=True,
        help="Exclude e2 carriers (default: True).",
    )

    return parser


def cmd_call(args: argparse.Namespace) -> None:
    """Execute the 'call' sub-command."""
    caller = APOECaller()

    fmt = args.format.lower()
    if fmt == "csv":
        results = caller.call_from_csv(args.input)
    elif fmt == "ped":
        results = caller.call_from_ped(args.input)
    elif fmt == "raw":
        results = caller.call_from_raw(args.input)
    else:
        print(f"Unsupported format: {fmt}", file=sys.stderr)
        sys.exit(1)

    df = caller.results_to_dataframe(results)

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Results written to {args.output} ({len(df):,} samples)")
    else:
        print(df.to_string(index=False))

    if args.summary:
        summary = caller.summarise(results)
        print(f"\n--- Summary ({summary.total_samples:,} samples) ---")
        for dip, count in sorted(summary.genotype_counts.items()):
            pct = (count / summary.total_samples * 100) if summary.total_samples else 0.0
            print(f"  {dip:12s}  {count:>6,}  ({pct:.1f}%)")
        if summary.indeterminate_count:
            print(f"  Indeterminate {summary.indeterminate_count:>6,}")


def cmd_feasibility(args: argparse.Namespace) -> None:
    """Execute the 'feasibility' sub-command."""
    estimator = APOEFeasibilityEstimator()
    report = estimator.estimate_from_csv(
        csv_path=args.input,
        target_diplotypes=args.targets,
        confidence_level=args.confidence,
    )
    print(estimator.format_report(report))


def cmd_stratify(args: argparse.Namespace) -> None:
    """Execute the 'stratify' sub-command."""
    import pandas as pd

    cohort = pd.read_csv(args.input)

    config = StratificationConfig(
        study_name=args.study,
        target_female_count=args.females,
        target_male_count=args.males,
        exclude_e2_carriers=args.exclude_e2,
    )

    stratifier = CohortStratifier()
    result = stratifier.stratify(cohort, config)
    files = stratifier.export_recall_lists(result, args.output_dir)

    print(f"Study: {config.study_name}")
    print(f"Total eligible: {result.total_eligible:,}")
    print(f"Excluded: {result.excluded_count:,}")
    print(f"Selected: {result.total_selected:,}")
    for f in files:
        print(f"  -> {f}")


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    dispatch = {
        "call": cmd_call,
        "feasibility": cmd_feasibility,
        "stratify": cmd_stratify,
    }

    handler = dispatch.get(args.command)
    if handler is None:
        parser.print_help()
        sys.exit(1)

    handler(args)


if __name__ == "__main__":
    main()
