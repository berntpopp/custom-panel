#!/bin/bash
# 
# Clean runner for custom-panel that suppresses all third-party deprecation warnings
# Usage: bash scripts/clean-run.sh [custom-panel arguments]
#
# This script provides a clean way to run custom-panel without seeing deprecation
# warnings from third-party dependencies like pypdf and cryptography.
#

export PYTHONWARNINGS="ignore::DeprecationWarning"
custom-panel "$@" 2>&1 | grep -vE "(CryptographyDeprecationWarning|ARC4 has been moved|from cryptography\.hazmat\.primitives\.ciphers\.algorithms import)"