#!/bin/bash
#SBATCH --job-name=apoe_genotyping
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=apoe_genotypes.%j.out
#SBATCH --error=apoe_genotypes.%j.err

# ==============================================================================
# Script Name: slurm_template.sh
# Description: Example SLURM submission script for running the APOE pipeline
#              on an HPC cluster.
# ==============================================================================

# 1. Load Modules (Adjust these names for your specific HPC environment)
module load plink/1.9
module load python/3.9

# 2. Configuration
# Path to the toolkit directory
TOOLKIT_DIR="$HOME/biobank-apoe-toolkit"
PIPELINE_SCRIPT="$TOOLKIT_DIR/run_apoe_pipeline.sh"

# Data Directories
DATA_DIR="/path/to/your/genotype/data"
OUTPUT_DIR="./results"
mkdir -p "$OUTPUT_DIR"

# 3. Run Pipeline for Multiple Batches
# Example: Loop through batch numbers 22, 25, 26
for BATCH in 22 25 26; do
    echo "Processing Batch $BATCH..."
    
    INPUT_FILE="$DATA_DIR/UKBBAffy_SAX_b${BATCH}"
    OUTPUT_PREFIX="$OUTPUT_DIR/batch${BATCH}_apoe"
    
    # Run the pipeline script
    bash "$PIPELINE_SCRIPT" "$INPUT_FILE" "$OUTPUT_PREFIX"
    
    echo "Finished Batch $BATCH"
done

echo "All jobs completed."
