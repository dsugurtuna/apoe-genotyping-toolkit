#!/bin/bash

# ==============================================================================
# Script Name: run_apoe_pipeline.sh
# Description: Automates the extraction of APOE SNPs from PLINK datasets and 
#              calls genotypes using the Python helper script.
# Usage:       ./run_apoe_pipeline.sh <input_bfile_prefix> <output_prefix>
# ==============================================================================

# --- Configuration ---
# Path to the Python script (assumes it's in the same directory as this script)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PYTHON_SCRIPT="$SCRIPT_DIR/apoe_caller.py"
SNP_LIST="$SCRIPT_DIR/snp_list.txt"

# Check for PLINK installation
if ! command -v plink &> /dev/null; then
    echo "Error: 'plink' command not found. Please load the module or install PLINK 1.9."
    exit 1
fi

# Check Arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_bfile_prefix> <output_prefix>"
    echo "Example: $0 /data/genotypes/batch28_vcf my_results/batch28_apoe"
    exit 1
fi

INPUT_BFILE=$1
OUTPUT_PREFIX=$2

echo "========================================================"
echo "APOE Genotyping Pipeline"
echo "Input:  $INPUT_BFILE"
echo "Output: $OUTPUT_PREFIX"
echo "========================================================"

# 1. Extract APOE SNPs (rs429358, rs7412)
echo "[Step 1] Extracting SNPs..."
plink \
  --bfile "$INPUT_BFILE" \
  --extract "$SNP_LIST" \
  --make-bed \
  --out "${OUTPUT_PREFIX}_temp" \
  --allow-extra-chr \
  --silent

if [ $? -ne 0 ]; then
    echo "Error: PLINK extraction failed."
    exit 1
fi

# 2. Convert to Compound Genotypes (.ped format)
# This ensures alleles are grouped (e.g., 'CT') for easier parsing
echo "[Step 2] Converting to compound genotypes..."
plink \
  --bfile "${OUTPUT_PREFIX}_temp" \
  --recode compound-genotypes \
  --out "${OUTPUT_PREFIX}_temp" \
  --silent

if [ $? -ne 0 ]; then
    echo "Error: PLINK recoding failed."
    exit 1
fi

# 3. Run Python Caller
echo "[Step 3] Calling Genotypes..."
python3 "$PYTHON_SCRIPT" \
  --input "${OUTPUT_PREFIX}_temp.ped" \
  --output "$OUTPUT_PREFIX"

# 4. Cleanup Temporary Files
echo "[Step 4] Cleaning up..."
rm "${OUTPUT_PREFIX}_temp"*

echo "========================================================"
echo "Pipeline Complete!"
echo "Results:"
echo "  - ${OUTPUT_PREFIX}.APOE_GENOTYPES.csv"
echo "  - ${OUTPUT_PREFIX}.APOE_SUMMARY.csv"
echo "========================================================"
