#!/bin/bash

# ==============================================================================
# Script Name: run_apoe_pipeline.sh
# Author:      Ugur Tuna
# Context:     Developed during tenure at NIHR BioResource (Cambridge).
# Disclaimer:  This script is a sanitized version for educational/portfolio use.
#              It contains no real patient data or internal infrastructure paths.
#
# Description: Automates the extraction of APOE SNPs from PLINK datasets and 
#              calls genotypes using the Python helper script.
# Usage:       ./run_apoe_pipeline.sh <input_bfile_prefix> <output_prefix>
# ==============================================================================

set -e  # Exit immediately if a command exits with a non-zero status

# --- Configuration ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PYTHON_SCRIPT="$SCRIPT_DIR/apoe_caller.py"
SNP_LIST="$SCRIPT_DIR/snp_list.txt"

# --- Helper Functions ---

print_usage() {
    echo "Usage: $0 <input_bfile_prefix> <output_prefix>"
    echo "Example: $0 /data/genotypes/batch28_vcf my_results/batch28_apoe"
    echo ""
    echo "Arguments:"
    echo "  input_bfile_prefix   Path to PLINK binary files (without .bed/.bim/.fam extension)"
    echo "  output_prefix        Path/Prefix for output files"
}

log_info() {
    echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

check_dependencies() {
    if ! command -v plink &> /dev/null; then
        log_error "'plink' command not found. Please load the module or install PLINK 1.9."
        exit 1
    fi
    if ! command -v python3 &> /dev/null; then
        log_error "'python3' command not found."
        exit 1
    fi
}

# --- Main Execution ---

# Check Arguments
if [ "$#" -ne 2 ]; then
    print_usage
    exit 1
fi

INPUT_BFILE=$1
OUTPUT_PREFIX=$2

check_dependencies

echo "========================================================"
echo "       APOE Genotyping Pipeline v1.0.0"
echo "========================================================"
log_info "Input:  $INPUT_BFILE"
log_info "Output: $OUTPUT_PREFIX"

# 1. Extract APOE SNPs (rs429358, rs7412)
log_info "Step 1: Extracting SNPs with PLINK..."
plink \
  --bfile "$INPUT_BFILE" \
  --extract "$SNP_LIST" \
  --make-bed \
  --out "${OUTPUT_PREFIX}_temp" \
  --allow-extra-chr \
  --silent

# 2. Convert to Compound Genotypes (.ped format)
log_info "Step 2: Converting to compound genotypes..."
plink \
  --bfile "${OUTPUT_PREFIX}_temp" \
  --recode compound-genotypes \
  --out "${OUTPUT_PREFIX}_temp" \
  --silent

# 3. Run Python Caller
log_info "Step 3: Calling Genotypes with Python..."
python3 "$PYTHON_SCRIPT" \
  --input "${OUTPUT_PREFIX}_temp.ped" \
  --output "$OUTPUT_PREFIX"

# 4. Cleanup Temporary Files
log_info "Step 4: Cleaning up temporary files..."
rm -f "${OUTPUT_PREFIX}_temp"*

echo "========================================================"
log_info "Pipeline Complete Successfully!"
echo "Results:"
echo "  - ${OUTPUT_PREFIX}.APOE_GENOTYPES.csv"
echo "  - ${OUTPUT_PREFIX}.APOE_SUMMARY.csv"
echo "========================================================"
