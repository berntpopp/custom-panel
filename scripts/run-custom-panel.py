#!/usr/bin/env python3
"""
Wrapper script to run custom-panel with deprecation warnings suppressed.

This script suppresses known deprecation warnings from dependencies
that are outside our control.
"""

import os
import warnings

# Set environment variable to suppress warnings at the interpreter level
os.environ["PYTHONWARNINGS"] = "ignore::DeprecationWarning"

# Also suppress ALL deprecation warnings globally at the earliest possible point
warnings.filterwarnings("ignore", category=DeprecationWarning)

# All other imports must come after warning suppression  # noqa: E402
from custom_panel.main import app  # noqa: E402

if __name__ == "__main__":
    app()
