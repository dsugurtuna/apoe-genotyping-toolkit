# Biobank APOE Genotyping Toolkit 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![PLINK](https://img.shields.io/badge/Tool-PLINK%201.9-red)](https://www.cog-genomics.org/plink/)

**Automated APOE Genotype Calling Pipeline for Biobank Scale Data.**

The **Biobank APOE Genotyping Toolkit** is a robust, HPC-ready pipeline designed to extract, process, and classify *Apolipoprotein E* (APOE) genotypes (e2/e2, e3/e4, etc.) from large-scale PLINK datasets. 

It standardizes the complex logic of mapping `rs429358` and `rs7412` haplotypes into clinically relevant APOE status, generating clean, audit-ready CSV reports.

---

## 🚀 Key Features

*   **⚡ Automated Extraction**: Seamlessly pulls relevant SNPs from massive PLINK binary files (`.bed`/`.bim`/`.fam`).
*   **🧬 Accurate Calling**: Implements standard haplotype mapping logic to determine all 6 APOE genotypes (e2/e2, e2/e3, e2/e4, e3/e3, e3/e4, e4/e4).
*   **📊 Summary Statistics**: Automatically calculates genotype frequencies and phenotype stratifications (Cases vs. Controls).
*   **🖥️ HPC Integration**: Includes SLURM templates for easy deployment on research clusters.

## 📂 Toolkit Contents

| File | Description |
| :--- | :--- |
| `run_apoe_pipeline.sh` | **Main Driver.** Orchestrates PLINK extraction and Python calling. |
| `apoe_caller.py` | **The Logic.** Python script that maps haplotypes to genotypes. |
| `merge_batches.sh` | **Utility.** Merges results from multiple batches into a master dataset. |
| `slurm_template.sh` | **HPC Job.** Example submission script for clusters. |

## 🛠️ Quick Start

### Prerequisites
*   **PLINK 1.9**: Must be installed or available via `module load plink`.
*   **Python 3.8+**: With `pandas` and `numpy` installed.

### Installation
```bash
git clone https://github.com/YOUR_ORG/biobank-apoe-toolkit.git
cd biobank-apoe-toolkit
chmod +x *.sh
```

### Usage

**1. Single Batch Run**
Run the pipeline directly on a PLINK dataset (do not include file extensions like .bed):
```bash
./run_apoe_pipeline.sh /path/to/data/batch28_genotypes ./output/batch28_results
```

**2. HPC Batch Processing**
Edit `slurm_template.sh` to point to your data directories, then submit:
```bash
sbatch slurm_template.sh
```

**3. Merging Results**
Combine outputs from multiple runs:
```bash
./merge_batches.sh final_merged_apoe.csv ./output/batch*_results.APOE_GENOTYPES.csv
```

## 🧬 Genotype Mapping Logic

The toolkit uses the standard mapping for the two SNPs:

| rs429358 | rs7412 | APOE Genotype |
| :--- | :--- | :--- |
| T | T | **e2/e2** |
| T | C | **e3/e3** |
| C | C | **e4/e4** |
| *Mixed* | *Mixed* | *Heterozygous combinations (e.g., e3/e4)* |

## 🤝 Contributing
We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details on how to submit pull requests.

## 📄 License
This project is licensed under the [MIT License](LICENSE).
