"""
Main entry point for the custom-panel CLI application.
"""

import os
import sys
import warnings
from typing import Any


# NUCLEAR OPTION: Completely suppress stderr warnings before any imports
class _StderrFilter:
    def __init__(self, original_stderr: Any) -> None:
        self.original_stderr = original_stderr

    def write(self, text: str) -> int:
        # Block specific deprecation warnings completely
        if any(
            pattern in text
            for pattern in [
                "CryptographyDeprecationWarning",
                "ARC4 has been moved",
                "from cryptography.hazmat.primitives.ciphers.algorithms import",
            ]
        ):
            return len(text)
        return self.original_stderr.write(text)

    def flush(self) -> None:
        self.original_stderr.flush()

    def __getattr__(self, name: str) -> Any:
        return getattr(self.original_stderr, name)


# Install stderr filter immediately
sys.stderr = _StderrFilter(sys.stderr)

# Set environment variable to suppress warnings at the interpreter level
os.environ.setdefault("PYTHONWARNINGS", "ignore::DeprecationWarning")

# Suppress ALL deprecation warnings globally at the earliest possible point
warnings.filterwarnings("ignore", category=DeprecationWarning)

from .cli import app  # noqa: E402

if __name__ == "__main__":
    app()
