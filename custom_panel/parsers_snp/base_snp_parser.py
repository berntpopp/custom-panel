"""
Base SNP parser abstract class.

This module provides the abstract base class for all SNP parsers in the custom-panel tool.
All SNP parsers should inherit from BaseSNPParser and implement the parse method.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import pandas as pd


class BaseSNPParser(ABC):
    """
    Abstract base class for SNP file parsers.

    This class defines the interface that all SNP parsers must implement.
    It provides common functionality and ensures consistent output format
    across all parser implementations.
    """

    def __init__(self, file_path: Path, config: dict[str, Any]):
        """
        Initialize the SNP parser.

        Args:
            file_path: Path to the SNP file to parse
            config: Configuration dictionary for the parser
        """
        self.file_path = Path(file_path)
        self.config = config
        self.name = config.get("name", "Unknown")

    @abstractmethod
    def parse(self) -> pd.DataFrame:
        """
        Parse the SNP file and return a standardized DataFrame.

        This method must be implemented by all concrete parser classes.
        The returned DataFrame must contain the following columns:
        - rsid: RS identifier (e.g., 'rs123456')
        - source: Source name/identifier
        - category: SNP category (e.g., 'identity', 'prs', 'manual')
        - additional metadata columns as needed

        Returns:
            Standardized DataFrame with parsed SNP data

        Raises:
            FileNotFoundError: If the file does not exist
            ValueError: If the file format is invalid
            Exception: For other parsing errors
        """

    def validate_file_exists(self) -> None:
        """
        Validate that the file exists and is readable.

        Raises:
            FileNotFoundError: If the file does not exist
        """
        if not self.file_path.exists():
            raise FileNotFoundError(f"SNP file not found: {self.file_path}")

        if not self.file_path.is_file():
            raise FileNotFoundError(f"Path is not a file: {self.file_path}")

    def validate_output_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Validate and standardize the output DataFrame format.

        Args:
            df: DataFrame to validate

        Returns:
            Validated and standardized DataFrame

        Raises:
            ValueError: If required columns are missing
        """
        # Check for required columns - either old format (rsid) or new format (snp)
        if "snp" in df.columns:
            # New format: snp is primary identifier, rsid is optional
            required_columns = ["snp", "source", "category"]
        else:
            # Old format: rsid is primary identifier
            required_columns = ["rsid", "source", "category"]

        # Check for required columns
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(
                f"Missing required columns in parser output: {missing_columns}",
            )

        df = df.copy()

        # Handle rsid column if present
        if "rsid" in df.columns:
            # Standardize rsid format (ensure they start with 'rs')
            df["rsid"] = df["rsid"].apply(self._standardize_rsid)

        # For new format, validate snp column instead of rsid
        if "snp" in df.columns:
            # Remove rows with empty snp identifiers
            df = df.dropna(subset=["snp"])
            df = df[df["snp"].str.strip() != ""]
        elif "rsid" in df.columns:
            # Old format: remove rows with empty rsids
            df = df.dropna(subset=["rsid"])
            df = df[df["rsid"].str.strip() != ""]

        # Ensure source and category are strings
        df["source"] = df["source"].astype(str)
        df["category"] = df["category"].astype(str)

        return df

    def _standardize_rsid(self, rsid: str | None) -> str | None:
        """
        Standardize rsID format to ensure it starts with 'rs'.

        Args:
            rsid: Raw rsID string

        Returns:
            Standardized rsID or None if invalid
        """
        if not rsid or pd.isna(rsid):
            return None

        rsid = str(rsid).strip()

        # Handle empty strings
        if not rsid:
            return None

        # If it already starts with 'rs', return as-is
        if rsid.lower().startswith("rs"):
            return rsid

        # If it's just numbers, add 'rs' prefix
        if rsid.isdigit():
            return f"rs{rsid}"

        # For other formats, return as-is but log a warning
        # (Could be chr:pos format or other identifier)
        return rsid

    def get_parser_info(self) -> dict[str, Any]:
        """
        Get information about this parser instance.

        Returns:
            Dictionary with parser metadata
        """
        return {
            "parser_class": self.__class__.__name__,
            "parser_name": self.name,
            "file_path": str(self.file_path),
            "file_exists": self.file_path.exists(),
            "config": self.config,
        }
