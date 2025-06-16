"""
Custom Panel: A modern Python tool for gene panel curation.

This package provides tools for aggregating gene information from multiple sources
and creating custom gene panels for clinical genomics applications.
"""

# AGGRESSIVE WARNING SUPPRESSION - MUST BE FIRST
import os
import sys
import warnings

# Set environment variable to suppress warnings at the interpreter level
os.environ.setdefault("PYTHONWARNINGS", "ignore::DeprecationWarning")

# Suppress ALL deprecation warnings globally and immediately
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Redirect stderr temporarily to suppress import-time warnings
class _WarningFilter:
    def __init__(self):
        self._original_stderr = sys.stderr
        self._buffer = []
        
    def write(self, text):
        if "CryptographyDeprecationWarning" in text or "ARC4 has been moved" in text:
            return  # Suppress these specific warnings
        self._original_stderr.write(text)
        
    def flush(self):
        self._original_stderr.flush()

# Temporarily replace stderr during package initialization
_temp_stderr = _WarningFilter()
sys.stderr = _temp_stderr

# Import everything else now
try:
    pass  # Package-level imports would go here
finally:
    # Restore original stderr
    sys.stderr = _temp_stderr._original_stderr

__version__ = "0.1.0"
