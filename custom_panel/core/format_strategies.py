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

    @abstractmethod
    def get_extension(self) -> str:
        """Get file extension for this format."""

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
            # Clean DataFrame for Parquet compatibility
            df_clean = self._clean_for_parquet(df)
            df_clean.to_parquet(path, index=index, engine=engine)
            logger.debug(f"Saved {len(df_clean)} records to Parquet: {path}")
        except Exception as e:
            logger.error(f"Failed to save Parquet file {path}: {e}")
            raise

    def _clean_for_parquet(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Clean DataFrame for Parquet compatibility.

        Handles mixed data types that cause Parquet conversion errors.
        """
        df_clean = df.copy()

        # Fix hg38_strand columns (mixed str/int types)
        for strand_col in ["hg38_strand"]:
            if strand_col in df_clean.columns:
                # Convert all values to string first, then standardize
                df_clean[strand_col] = df_clean[strand_col].astype(str)
                # Replace empty strings and 'nan' with None for consistency
                df_clean[strand_col] = df_clean[strand_col].replace(["", "nan"], None)
                # Convert numeric strings back to integers where appropriate
                df_clean[strand_col] = df_clean[strand_col].apply(
                    lambda x: int(x) if x and x.isdigit() else x,
                )

        return df_clean

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
        snp_data = kwargs.get("snp_data")
        regions_data = kwargs.get("regions_data")

        try:
            # Check if this is a multi-sheet Excel file with SNP or regions data
            if (snp_data and isinstance(snp_data, dict)) or (
                regions_data and isinstance(regions_data, dict)
            ):
                self._save_multi_sheet_excel(
                    df, path, snp_data, regions_data, index, engine,
                )
            else:
                # Single sheet Excel
                df.to_excel(path, index=index, engine=engine, sheet_name=sheet_name)
                logger.debug(f"Saved {len(df)} records to Excel: {path}")
        except Exception as e:
            logger.error(f"Failed to save Excel file {path}: {e}")
            raise

    def _save_multi_sheet_excel(
        self,
        gene_df: pd.DataFrame,
        path: Path,
        snp_data: dict[str, pd.DataFrame] | None,
        regions_data: dict[str, pd.DataFrame] | None,
        index: bool,
        engine: Any,
    ) -> None:
        """Save multi-sheet Excel with genes, SNPs, and regions."""
        from rich.console import Console

        console = Console()

        with pd.ExcelWriter(path, engine=engine) as writer:
            # Save gene panel as first sheet
            gene_df.to_excel(writer, sheet_name="Gene_Panel", index=index)
            logger.debug(f"Saved {len(gene_df)} genes to Gene_Panel sheet")

            # Process SNP data if available
            if snp_data and isinstance(snp_data, dict):
                # Combine all SNPs for master SNP sheet
                all_snps = []
                for _, snp_df in snp_data.items():
                    if not snp_df.empty:
                        snp_df_copy = snp_df.copy()
                        # Note: snp_type info already preserved in category column
                        all_snps.append(snp_df_copy)

                if all_snps:
                    master_snp_df = pd.concat(all_snps, ignore_index=True)
                    master_snp_df.to_excel(writer, sheet_name="All_SNPs", index=index)
                    logger.debug(f"Saved {len(master_snp_df)} SNPs to All_SNPs sheet")
                    console.print(
                        f"[blue]Added All_SNPs sheet with {len(master_snp_df)} SNPs[/blue]",
                    )

                    # Add individual SNP type sheets
                    for snp_type, snp_df in snp_data.items():
                        if not snp_df.empty:
                            # Create valid sheet name (Excel sheet names <= 31 chars)
                            sheet_name = f"SNPs_{snp_type}"[:31]
                            snp_df.to_excel(writer, sheet_name=sheet_name, index=index)
                            logger.debug(
                                f"Saved {len(snp_df)} SNPs to {sheet_name} sheet",
                            )
                            console.print(
                                f"[blue]Added {sheet_name} sheet with "
                                f"{len(snp_df)} SNPs[/blue]",
                            )

            # Process regions data if available
            if regions_data and isinstance(regions_data, dict):
                # Combine all regions for master regions sheet
                all_regions = []
                for region_type, region_df in regions_data.items():
                    if not region_df.empty:
                        region_df_copy = region_df.copy()
                        region_df_copy["region_type"] = region_type
                        all_regions.append(region_df_copy)

                if all_regions:
                    master_regions_df = pd.concat(all_regions, ignore_index=True)
                    master_regions_df.to_excel(
                        writer, sheet_name="All_Regions", index=index,
                    )
                    logger.debug(
                        f"Saved {len(master_regions_df)} regions to All_Regions sheet",
                    )
                    console.print(
                        f"[blue]Added All_Regions sheet with {len(master_regions_df)} regions[/blue]",
                    )

                    # Add individual region type sheets
                    for region_type, region_df in regions_data.items():
                        if not region_df.empty:
                            # Create valid sheet name (Excel sheet names <= 31 chars)
                            sheet_name = f"Regions_{region_type}"[:31]
                            region_df.to_excel(
                                writer, sheet_name=sheet_name, index=index,
                            )
                            logger.debug(
                                f"Saved {len(region_df)} regions to {sheet_name} sheet",
                            )
                            console.print(
                                f"[blue]Added {sheet_name} sheet with "
                                f"{len(region_df)} regions[/blue]",
                            )

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


class BedStrategy(FormatStrategy):
    """Strategy for saving DataFrames as BED files."""

    def save(self, df: pd.DataFrame, path: Path, **kwargs: Any) -> None:
        """Save DataFrame as BED file."""
        path.parent.mkdir(parents=True, exist_ok=True)

        # BED format expects tab-separated values without headers
        # Standard BED format: chrom, start, end, name, score, strand
        try:
            df.to_csv(path, sep="\t", header=False, index=False)
            logger.debug(f"Saved {len(df)} records to BED: {path}")
        except Exception as e:
            logger.error(f"Failed to save BED file {path}: {e}")
            raise

    def get_extension(self) -> str:
        return "bed"


class FormatStrategyFactory:
    """Factory for creating format strategy instances."""

    _strategies: dict[str, type[FormatStrategy]] = {
        "parquet": ParquetStrategy,
        "csv": CSVStrategy,
        "excel": ExcelStrategy,
        "xlsx": ExcelStrategy,  # Alias for excel
        "json": JSONStrategy,
        "bed": BedStrategy,
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
                f"Unsupported format: {format_name}. Available: {available}",
            )

        return cls._strategies[format_name]()

    @classmethod
    def get_supported_formats(cls) -> list[str]:
        """Get list of supported format names."""
        return list(cls._strategies.keys())

    @classmethod
    def register_strategy(
        cls, format_name: str, strategy_class: type[FormatStrategy],
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
        self, df: pd.DataFrame, path: Path, format_name: str, **kwargs: Any,
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
        self, base_path: Path, filename_base: str, format_name: str,
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
    df: pd.DataFrame, path: Path, format_name: str, **kwargs: Any,
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
