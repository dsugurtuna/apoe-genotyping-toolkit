#!/usr/bin/env python3
"""
APOE Genotype Caller
====================

A robust bioinformatics utility for determining Apolipoprotein E (APOE) genotypes 
from phased or unphased PLINK data.

This script implements the standard clinical mapping logic for rs429358 and rs7412 
variants to classify samples into e2/e2, e3/e3, e4/e4, and heterozygous combinations.

Features:
- Handles PLINK .ped format input
- Generates detailed summary statistics
- Stratifies by phenotype (Case/Control)
- Type-safe and error-tolerant

Author: Ugur Tuna
License: MIT
"""

import argparse
import sys
from typing import List, Dict, Optional
from functools import reduce
import pandas as pd

# ==============================================================================
# Configuration & Constants
# ==============================================================================

# Standard PLINK .ped columns
PED_HEADER: List[str] = ["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]

# The specific SNPs required for APOE classification
REQUIRED_SNPS: List[str] = ["rs429358", "rs7412"]

# Mapping logic: (rs429358_allele + rs7412_allele) -> Genotype
# Note: Alleles must be consistent with the reference genome used (usually forward strand)
APOE_MAPPING: Dict[str, str] = {
    # Homozygous
    'TT_TT': 'e2/e2',
    'TT_CC': 'e3/e3',
    'CC_CC': 'e4/e4',
    
    # Heterozygous
    'CT_TT': 'e1/e2', 'TC_TT': 'e1/e2', # Rare e1 variants
    'CC_TT': 'e1/e1',
    
    'TT_TC': 'e2/e3', 'TT_CT': 'e2/e3',
    'TC_CC': 'e3/e4', 'CT_CC': 'e3/e4',
    
    # Ambiguous / Complex Haplotypes (e2/e4 vs e1/e3)
    # Often clinically reported as e2/e4 due to higher frequency
    'TC_TC': 'e2/e4', 'CT_CT': 'e2/e4', 
    'TC_CT': 'e2/e4', 'CT_TC': 'e2/e4',
    
    'CC_CT': 'e1/e4', 'CC_TC': 'e1/e4'
}

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Determine APOE genotypes from PLINK .ped data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", "-i", required=True, help="Input .ped file path")
    parser.add_argument("--output", "-o", required=True, help="Output prefix for CSV files")
    return parser.parse_args()

def load_ped_file(filepath: str) -> pd.DataFrame:
    """
    Load a PLINK .ped file into a pandas DataFrame.
    
    Args:
        filepath: Path to the .ped file.
        
    Returns:
        pd.DataFrame: Loaded data with appropriate column names.
    """
    # The file is expected to have standard PED columns + the 2 extracted SNPs
    header = PED_HEADER + REQUIRED_SNPS
    
    try:
        # Use 'sep' for whitespace handling (PLINK default)
        df = pd.read_csv(filepath, sep=r'\s+', names=header, header=None, engine='python')
        return df
    except FileNotFoundError:
        print(f"Error: Input file '{filepath}' not found.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Input file '{filepath}' is empty.")
        sys.exit(1)

def determine_genotypes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Map SNP combinations to APOE genotypes.
    """
    # Create a combined key for mapping
    # We cast to string to handle potential numeric allele coding (though usually A/C/G/T)
    df['haplotype_key'] = (
        df['rs429358'].astype(str) + "_" + df['rs7412'].astype(str)
    )

    # Map and fill unknowns
    df['APOE_GENOTYPE'] = df['haplotype_key'].map(APOE_MAPPING).fillna('Indeterminate')
    
    return df

def generate_summary_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics stratified by phenotype.
    """
    # Helper to calculate counts and percentages
    def get_counts(sub_df: pd.DataFrame, label: str) -> pd.DataFrame:
        if sub_df.empty:
            return pd.DataFrame(columns=['APOE_GENOTYPE', f'{label}_COUNT', f'{label}_PERCENT'])
            
        c = sub_df['APOE_GENOTYPE'].value_counts().reset_index()
        c.columns = ['APOE_GENOTYPE', f'{label}_COUNT']
        c[f'{label}_PERCENT'] = (c[f'{label}_COUNT'] / len(sub_df) * 100).round(2)
        return c

    # 1. Overall Counts
    total_counts = get_counts(df, 'TOTAL')

    # 2. Stratify by Phenotype (PLINK standard: -9=Missing, 1=Control, 2=Case)
    missing_pheno = df[df['PHENO'] == -9]
    controls = df[df['PHENO'] == 1]
    cases = df[df['PHENO'] == 2]

    m_counts = get_counts(missing_pheno, 'MISSING_PHENO')
    ctrl_counts = get_counts(controls, 'CONTROLS')
    case_counts = get_counts(cases, 'CASES')

    # 3. Merge all summaries
    dfs = [total_counts, m_counts, ctrl_counts, case_counts]
    summary_df = reduce(
        lambda left, right: pd.merge(left, right, on='APOE_GENOTYPE', how='outer'), 
        dfs
    )
    
    return summary_df.fillna(0)

def main():
    args = parse_args()
    print(f"--- APOE Genotype Caller Started ---")
    print(f"Input: {args.input}")

    # 1. Load Data
    df = load_ped_file(args.input)
    print(f"Loaded {len(df)} samples.")

    # 2. Call Genotypes
    df = determine_genotypes(df)
    
    # 3. Generate Statistics
    summary_df = generate_summary_stats(df)

    # 4. Save Outputs
    # Drop intermediate columns for clean output
    output_cols = PED_HEADER + ['APOE_GENOTYPE']
    final_df = df[output_cols]
    
    geno_out = f"{args.output}.APOE_GENOTYPES.csv"
    summary_out = f"{args.output}.APOE_SUMMARY.csv"
    
    final_df.to_csv(geno_out, index=False)
    summary_df.to_csv(summary_out, index=False)

    print(f"[SUCCESS] Genotypes saved to: {geno_out}")
    print(f"[SUCCESS] Summary saved to:    {summary_out}")

if __name__ == "__main__":
    main()
