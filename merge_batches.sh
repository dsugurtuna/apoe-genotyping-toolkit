#!/bin/bash

# ==============================================================================
# Script Name: merge_batches.sh
# Description: Merges multiple APOE genotype CSV files into a single master file.
# Usage:       ./merge_batches.sh <output_filename.csv> <input_file1.csv> [input_file2.csv ...]
# ==============================================================================

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <output_filename.csv> <input_file1.csv> [input_file2.csv ...]"
    exit 1
fi

OUTPUT_FILE=$1
shift # Remove the first argument (output file) from the list

# Create the header from the first file
head -n 1 "$1" > "$OUTPUT_FILE"

echo "Merging files into $OUTPUT_FILE..."

# Loop through all remaining arguments (input files)
for INPUT_FILE in "$@"; do
    if [ -f "$INPUT_FILE" ]; then
        echo "  Appending: $INPUT_FILE"
        # Append content, skipping the header (line 1)
        tail -n +2 "$INPUT_FILE" >> "$OUTPUT_FILE"
    else
        echo "  Warning: File '$INPUT_FILE' not found. Skipping."
    fi
done

echo "-------------------------------------------------------"
echo "Merge Complete."
echo "Total lines (including header): $(wc -l < "$OUTPUT_FILE")"
