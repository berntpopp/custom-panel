"""
Main entry point for the custom-panel CLI application.
"""

import warnings

# Suppress known deprecation warnings before any other imports
warnings.filterwarnings(
    "ignore", message=".*ARC4 has been moved.*", category=DeprecationWarning
)
warnings.filterwarnings(
    "ignore", message=".*'BaseCommand' is deprecated.*", category=DeprecationWarning
)

from .cli import app

if __name__ == "__main__":
    app()
