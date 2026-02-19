"""
APOE Genotyping Toolkit
=======================

A production-grade bioinformatics toolkit for determining Apolipoprotein E (APOE)
genotypes from large-scale genomic datasets. Supports genotype calling, feasibility
estimation, and stratified recall list generation for clinical studies.

Author: Ugur Tuna
"""

__version__ = "2.0.0"

from apoe_toolkit.caller import APOECaller
from apoe_toolkit.feasibility import APOEFeasibilityEstimator
from apoe_toolkit.stratifier import CohortStratifier

__all__ = ["APOECaller", "APOEFeasibilityEstimator", "CohortStratifier"]
