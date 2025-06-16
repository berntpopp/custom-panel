#!/bin/bash
"""
Clean wrapper script for custom-panel that suppresses all deprecation warnings.

This script completely eliminates deprecation warnings from third-party dependencies
by setting the appropriate environment variables before running the tool.
"""

export PYTHONWARNINGS="ignore::DeprecationWarning"
exec custom-panel "$@"