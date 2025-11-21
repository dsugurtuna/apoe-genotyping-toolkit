# Biobank APOE Genotyping Toolkit 🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![PLINK](https://img.shields.io/badge/Tool-PLINK%201.9-red)](https://www.cog-genomics.org/plink/)
[![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics-green.svg)]()

**A High-Performance, Automated Pipeline for Clinical APOE Genotyping at Biobank Scale.**

The **Biobank APOE Genotyping Toolkit** is a production-grade bioinformatics pipeline designed to standardize the extraction and classification of *Apolipoprotein E* (APOE) genotypes from large-scale genomic datasets (PLINK).

It is engineered for **reproducibility**, **scalability**, and **clinical accuracy**, making it an essential tool for Alzheimer's research and population genetics.

---

## 🌟 Why This Matters

The APOE gene is the strongest genetic risk factor for late-onset Alzheimer's disease. Accurately determining a participant's genotype (e.g., `e3/e4` vs `e3/e3`) is critical for:
*   **Clinical Trials**: Stratifying patients by risk profile.
*   **GWAS Studies**: Adjusting for population structure and risk covariates.
*   **Precision Medicine**: Tailoring interventions based on genetic susceptibility.

This toolkit solves the challenge of manually parsing complex haplotype combinations (`rs429358` + `rs7412`) across millions of samples.

## 🚀 Key Features

*   **⚡ HPC-Optimized**: Built to run on Slurm/SGE clusters, processing thousands of samples in minutes.
*   **🛡️ Robust Classification**: Implements standard clinical mapping logic to resolve all 6 genotype combinations, including rare variants.
*   **📊 Automated Reporting**: Generates audit-ready CSVs with phenotype stratification (Cases vs. Controls).
*   **🔧 Modular Design**: Decoupled extraction (Bash/PLINK) and logic (Python) layers for easy maintenance.

## 🧬 Genotype Mapping Logic

The toolkit implements the following standard haplotype mapping:

| rs429358 (112) | rs7412 (158) | APOE Genotype | Risk Profile |
| :--- | :--- | :--- | :--- |
| **T** | **T** | **e2/e2** | Lower Risk |
| **T** | **C** | **e3/e3** | Neutral (Benchmark) |
| **C** | **C** | **e4/e4** | High Risk (Alzheimer's) |
| *Heterozygous* | *Heterozygous* | *e.g., e3/e4* | Mixed Risk |

## 📂 Repository Structure

```text
.
├── run_apoe_pipeline.sh   # 🚀 Main Driver: Orchestrates PLINK & Python
├── apoe_caller.py         # 🧠 Core Logic: Genotype classification engine
├── merge_batches.sh       # 🔄 Utility: Merges multi-batch outputs
├── slurm_template.sh      # ⚡ HPC: Ready-to-use cluster submission script
├── snp_list.txt           # ⚙️ Config: Target SNP definitions
└── requirements.txt       # 📦 Dependencies
```

## 🛠️ Quick Start

### Prerequisites
*   **PLINK 1.9**: [Download](https://www.cog-genomics.org/plink/)
*   **Python 3.8+**: `pip install -r requirements.txt`

### Installation
```bash
git clone https://github.com/dsugurtuna/biobank-apoe-toolkit.git
cd biobank-apoe-toolkit
chmod +x *.sh
```

### Usage

**1. Standard Execution**
Run the pipeline on a local PLINK dataset:
```bash
./run_apoe_pipeline.sh /path/to/data/batch28_genotypes ./results/batch28
```

**2. HPC Batch Processing (Slurm)**
Submit a job to process multiple batches in parallel:
```bash
sbatch slurm_template.sh
```

**3. Merging Results**
Combine outputs from multiple cohorts into a master dataset:
```bash
./merge_batches.sh master_apoe_calls.csv ./results/*.APOE_GENOTYPES.csv
```

## 🤝 Contributing
Contributions are welcome! We follow strict coding standards to ensure clinical accuracy. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## 📄 License
This project is licensed under the [MIT License](LICENSE).

---
*Designed for the scientific community to advance Alzheimer's research through transparent and reproducible data science.*
