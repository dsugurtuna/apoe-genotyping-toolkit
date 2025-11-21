#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
from functools import reduce

# ==============================================================================
# Script Name: apoe_caller.py
# Description: Determines APOE genotypes (e.g., e3/e4) from PLINK .ped files
#              containing rs429358 and rs7412.
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Determine APOE genotypes from PLINK .ped data.")
    parser.add_argument("--input", "-i", required=True, help="Input .ped file path")
    parser.add_argument("--output", "-o", required=True, help="Output prefix for CSV files")
    return parser.parse_args()

def main():
    args = parse_args()

    # PLINK .ped standard columns + the 2 SNPs we extracted
    # Note: Ensure your PLINK extraction order matches this, or use --recode compound-genotypes
    header_text = ["FID","IID","PAT","MAT","SEX","PHENO","rs429358","rs7412"]
    
    try:
        df = pd.read_csv(args.input, sep=" ", names=header_text, header=None)
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found.")
        exit(1)

    # Create a combined column for mapping
    df['rs429358_rs7412'] = df['rs429358'].astype(str) + "_" + df['rs7412'].astype(str)

    # APOE Genotype Mapping Dictionary
    # Based on rs429358 (C/T) and rs7412 (C/T) alleles
    apoe_dict = {
        'CC_TT': 'e1/e1', 'CT_TT': 'e1/e2','TC_TT': 'e1/e2',
        'CC_CT': 'e1/e4','CC_TC': 'e1/e4',
        'TT_TT': 'e2/e2',
        'TT_TC': 'e2/e3','TT_CT': 'e2/e3',
        'TC_TC': 'e2/e4 or e1/e3','CT_CT': 'e2/e4 or e1/e3','TC_CT': 'e2/e4 or e1/e3','CT_TC': 'e2/e4 or e1/e3',
        'TT_CC': 'e3/e3',
        'TC_CC': 'e3/e4','CT_CC': 'e3/e4',
        'CC_CC': 'e4/e4'
    }

    df['APOE_GENOTYPE'] = df['rs429358_rs7412'].map(apoe_dict).fillna('unknown')

    # Prepare Genotype Output
    subset = df.drop(columns=['PAT','MAT','rs429358','rs7412'])
    
    # --- Generate Summary Statistics ---
    counts = subset['APOE_GENOTYPE'].value_counts().reset_index()
    counts.columns = ['APOE_GENOTYPE','TOTAL_COUNT']
    counts['TOTAL_PERCENT'] = counts['TOTAL_COUNT']/subset.shape[0]*100

    # Subsets based on Phenotype (if available in .ped)
    # -9 = Missing, 1 = Control, 2 = Case
    missing   = subset[subset['PHENO']==-9]
    controls  = subset[subset['PHENO']== 1]
    cases     = subset[subset['PHENO']== 2]

    def get_counts(df, label):
        c = df['APOE_GENOTYPE'].value_counts().reset_index()
        c.columns = ['APOE_GENOTYPE', f'{label}_COUNT']
        c[f'{label}_PERCENT'] = c[f'{label}_COUNT']/df.shape[0]*100 if df.shape[0]>0 else 0
        return c

    m_counts  = get_counts(missing, 'MISSING_PHENO')
    ctrl_cts  = get_counts(controls,'CONTROLS')
    case_cts  = get_counts(cases,   'CASES')

    # Merge all summaries
    dataframes = [counts, m_counts, ctrl_cts, case_cts]
    summary_df = reduce(lambda left,right: pd.merge(left,right,on='APOE_GENOTYPE',how='left'), dataframes)

    # Save Outputs
    geno_out = args.output + ".APOE_GENOTYPES.csv"
    summary_out = args.output + ".APOE_SUMMARY.csv"
    
    subset.to_csv(geno_out, index=False)
    summary_df.to_csv(summary_out, index=False)

    print(f"[SUCCESS] Genotypes saved to: {geno_out}")
    print(f"[SUCCESS] Summary saved to:    {summary_out}")

if __name__ == "__main__":
    main()
