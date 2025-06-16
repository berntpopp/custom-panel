#!/usr/bin/env python3
"""
Wrapper script to run custom-panel with deprecation warnings suppressed.

This script suppresses known deprecation warnings from dependencies
that are outside our control.
"""

import sys
import warnings

# Suppress deprecation warnings before any other imports
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Now import and run the main application
from custom_panel.main import app

if __name__ == "__main__":
    app()
