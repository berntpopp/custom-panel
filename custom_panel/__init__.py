"""
Custom Panel: A modern Python tool for gene panel curation.

This package provides tools for aggregating gene information from multiple sources
and creating custom gene panels for clinical genomics applications.
"""

import warnings

# Suppress known deprecation warnings at package import time
warnings.filterwarnings(
    "ignore", message=".*ARC4 has been moved.*", category=DeprecationWarning
)
warnings.filterwarnings(
    "ignore", message=".*'BaseCommand' is deprecated.*", category=DeprecationWarning
)

__version__ = "0.1.0"
