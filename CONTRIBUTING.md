# Contributing to Biobank APOE Toolkit

We are excited that you are interested in contributing! This toolkit is designed to help researchers standardize APOE calling, and your improvements can help the entire community.

## How to Contribute

1.  **Fork the Repository**: Create your own copy of the project.
2.  **Create a Branch**: `git checkout -b feature/new-mapping-logic`
3.  **Make Changes**: Improve the Python logic or Bash scripts.
4.  **Test**: Ensure the pipeline still runs on standard PLINK files.
5.  **Submit a Pull Request**: Describe your changes and why they are needed.

## Development Guidelines

*   **Python**: Use `pandas` for data manipulation. Keep the code compatible with Python 3.8+.
*   **Bash**: Use clear variable names and comments.
*   **HPC Compatibility**: Avoid hardcoding paths specific to one cluster. Use variables.

## Reporting Issues

Found a bug in the genotype mapping? Please open an issue immediately with the specific SNP combination that failed.
