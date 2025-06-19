"""
File format strategies to eliminate repetition in data saving operations.

This module implements the Strategy pattern for different file formats,
providing a clean interface for saving DataFrames in various formats.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd
else:
    import pandas as pd

logger = logging.getLogger(__name__)


class FormatStrategy(ABC):
    """Abstract base class for file format strategies."""

    @abstractmethod
    def save(self, df: pd.DataFrame, path: Path, **kwargs: Any) -> None:
        """
        Save DataFrame to file in specific format.

        Args:
            df: DataFrame to save
            path: Output file path
            **kwargs: Additional format-specific options
        """
        pass

    @abstractmethod
    def get_extension(self) -> str:
        """Get file extension for this format."""
        pass

    def prepare_path(self, base_path: Path, filename: str) -> Path:
        """
        Prepare full file path with correct extension.

        Args:
            base_path: Base directory path
            filename: Base filename without extension

        Returns:
            Full path with correct extension
        """
        return base_path / f"{filename}.{self.get_extension()}"


class ParquetStrategy(FormatStrategy):
    """Strategy for saving DataFrames as Parquet files."""

    def save(self, df: pd.DataFrame, path: Path, **kwargs: Any) -> None:
        """Save DataFrame as Parquet file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        engine = kwargs.get("engine", "pyarrow")
        index = kwargs.get("index", False)

        try:
            df.to_parquet(path, index=index, engine=engine)
            logger.debug(f"Saved {len(df)} records to Parquet: {path}")
        except Exception as e:
            logger.error(f"Failed to save Parquet file {path}: {e}")
            raise

    def get_extension(self) -> str:
        return "parquet"


class CSVStrategy(FormatStrategy):
    """Strategy for saving DataFrames as CSV files."""

    def save(self, df: pd.DataFrame, path: Path, **kwargs: Any) -> None:
        """Save DataFrame as CSV file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        index = kwargs.get("index", False)
        encoding = kwargs.get("encoding", "utf-8")

        try:
            df.to_csv(path, index=index, encoding=encoding)
            logger.debug(f"Saved {len(df)} records to CSV: {path}")
        except Exception as e:
            logger.error(f"Failed to save CSV file {path}: {e}")
            raise

    def get_extension(self) -> str:
        return "csv"


class ExcelStrategy(FormatStrategy):
    """Strategy for saving DataFrames as Excel files."""

    def save(self, df: pd.DataFrame, path: Path, **kwargs: Any) -> None:
        """Save DataFrame as Excel file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        index = kwargs.get("index", False)
        engine = kwargs.get("engine", "openpyxl")
        sheet_name = kwargs.get("sheet_name", "Sheet1")

        try:
            df.to_excel(path, index=index, engine=engine, sheet_name=sheet_name)
            logger.debug(f"Saved {len(df)} records to Excel: {path}")
        except Exception as e:
            logger.error(f"Failed to save Excel file {path}: {e}")
            raise

    def get_extension(self) -> str:
        return "xlsx"


class JSONStrategy(FormatStrategy):
    """Strategy for saving DataFrames as JSON files."""

    def save(self, df: pd.DataFrame, path: Path, **kwargs: Any) -> None:
        """Save DataFrame as JSON file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        orient = kwargs.get("orient", "records")
        indent = kwargs.get("indent", 2)

        try:
            df.to_json(path, orient=orient, indent=indent)
            logger.debug(f"Saved {len(df)} records to JSON: {path}")
        except Exception as e:
            logger.error(f"Failed to save JSON file {path}: {e}")
            raise

    def get_extension(self) -> str:
        return "json"


class FormatStrategyFactory:
    """Factory for creating format strategy instances."""

    _strategies: dict[str, type[FormatStrategy]] = {
        "parquet": ParquetStrategy,
        "csv": CSVStrategy,
        "excel": ExcelStrategy,
        "xlsx": ExcelStrategy,  # Alias for excel
        "json": JSONStrategy,
    }

    @classmethod
    def create(cls, format_name: str) -> FormatStrategy:
        """
        Create a format strategy instance.

        Args:
            format_name: Name of the format (parquet, csv, excel, json)

        Returns:
            Format strategy instance

        Raises:
            ValueError: If format is not supported
        """
        format_name = format_name.lower()

        if format_name not in cls._strategies:
            available = ", ".join(cls._strategies.keys())
            raise ValueError(
                f"Unsupported format: {format_name}. Available: {available}"
            )

        return cls._strategies[format_name]()

    @classmethod
    def get_supported_formats(cls) -> list[str]:
        """Get list of supported format names."""
        return list(cls._strategies.keys())

    @classmethod
    def register_strategy(
        cls, format_name: str, strategy_class: type[FormatStrategy]
    ) -> None:
        """
        Register a new format strategy.

        Args:
            format_name: Name of the format
            strategy_class: Strategy class to register
        """
        cls._strategies[format_name.lower()] = strategy_class


class DataFrameSaver:
    """High-level interface for saving DataFrames in multiple formats."""

    def __init__(self) -> None:
        self.factory = FormatStrategyFactory()

    def save(
        self, df: pd.DataFrame, path: Path, format_name: str, **kwargs: Any
    ) -> None:
        """
        Save DataFrame using specified format strategy.

        Args:
            df: DataFrame to save
            path: Output file path
            format_name: Format to use (parquet, csv, excel, json)
            **kwargs: Additional format-specific options
        """
        strategy = self.factory.create(format_name)
        strategy.save(df, path, **kwargs)

    def save_multiple_formats(
        self,
        df: pd.DataFrame,
        base_path: Path,
        filename_base: str,
        formats: list[str],
        **kwargs: Any,
    ) -> dict[str, Path]:
        """
        Save DataFrame in multiple formats.

        Args:
            df: DataFrame to save
            base_path: Base directory for output files
            filename_base: Base filename without extension
            formats: List of formats to save in
            **kwargs: Additional format-specific options

        Returns:
            Dictionary mapping format names to saved file paths
        """
        saved_files = {}

        for format_name in formats:
            try:
                strategy = self.factory.create(format_name)
                file_path = strategy.prepare_path(base_path, filename_base)
                strategy.save(df, file_path, **kwargs)
                saved_files[format_name] = file_path
            except Exception as e:
                logger.error(f"Failed to save {format_name} format: {e}")
                # Continue with other formats

        return saved_files

    def get_file_path(
        self, base_path: Path, filename_base: str, format_name: str
    ) -> Path:
        """
        Get the file path for a specific format without saving.

        Args:
            base_path: Base directory
            filename_base: Base filename without extension
            format_name: Format name

        Returns:
            Full file path with correct extension
        """
        strategy = self.factory.create(format_name)
        return strategy.prepare_path(base_path, filename_base)


# Convenience function for backward compatibility
def save_dataframe(
    df: pd.DataFrame, path: Path, format_name: str, **kwargs: Any
) -> None:
    """
    Convenience function to save a DataFrame in the specified format.

    Args:
        df: DataFrame to save
        path: Output file path
        format_name: Format to use
        **kwargs: Additional format-specific options
    """
    saver = DataFrameSaver()
    saver.save(df, path, format_name, **kwargs)
