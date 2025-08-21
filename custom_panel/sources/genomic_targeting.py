"""
Genomic targeting flags data source.

This module processes genomic targeting flags from external files to mark genes
for complete genomic targeting. The flags are applied as a final enhancement step
in the pipeline after gene symbol standardization.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.config_manager import ConfigManager

logger = logging.getLogger(__name__)


def fetch_genomic_targeting_flags(config: dict[str, Any]) -> dict[str, bool]:
    """
    Fetch genomic targeting flags from configured file.

    Args:
        config: Configuration dictionary

    Returns:
        Dictionary mapping gene symbols to targeting flags (True/False)
    """
    config_manager = ConfigManager(config)
    targeting_config = config_manager.get_nested("genomic_targeting", default={})

    if not targeting_config.get("enabled", False):
        logger.info("Genomic targeting flags are disabled")
        return {}

    file_path = targeting_config.get("file_path")
    if not file_path:
        logger.warning("No genomic targeting file path configured")
        return {}

    file_path = Path(file_path)

    # Handle missing file gracefully
    allow_missing_file = targeting_config.get("allow_missing_file", True)
    if not file_path.exists():
        if allow_missing_file:
            logger.warning(
                f"Genomic targeting file not found: {file_path} - using defaults",
            )
            return {}
        else:
            raise FileNotFoundError(f"Genomic targeting file not found: {file_path}")

    logger.info(f"Loading genomic targeting flags from {file_path}")

    try:
        targeting_flags = process_targeting_file(file_path, targeting_config)
        logger.info(
            f"Successfully loaded {len(targeting_flags)} genomic targeting flags",
        )
        return targeting_flags
    except Exception as e:
        logger.error(f"Failed to process genomic targeting file: {e}")
        if allow_missing_file:
            logger.warning("Continuing with empty targeting flags due to error")
            return {}
        else:
            raise


def process_targeting_file(
    file_path: Path,
    targeting_config: dict[str, Any],
) -> dict[str, bool]:
    """
    Process genomic targeting file and return gene -> targeting flag mapping.

    Args:
        file_path: Path to the targeting file
        targeting_config: Configuration for targeting processing

    Returns:
        Dictionary mapping gene symbols to targeting flags
    """
    gene_column = targeting_config.get("gene_column", "gene_symbol")
    targeting_column = targeting_config.get("targeting_column", "targeting")
    default_value = targeting_config.get("default_value", False)

    # Read file based on extension
    file_extension = file_path.suffix.lower()

    try:
        if file_extension in [".xlsx", ".xls"]:
            df = pd.read_excel(file_path)
        elif file_extension == ".csv":
            df = pd.read_csv(file_path)
        elif file_extension in [".txt", ".tsv"]:
            df = pd.read_csv(file_path, sep="\t")
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")

        if df.empty:
            logger.warning("Genomic targeting file is empty")
            return {}

    except Exception as e:
        raise ValueError(f"Error reading genomic targeting file: {e}") from e

    # Validate required columns
    if gene_column not in df.columns:
        raise ValueError(f"Gene column '{gene_column}' not found in targeting file")

    if targeting_column not in df.columns:
        raise ValueError(
            f"Targeting column '{targeting_column}' not found in targeting file",
        )

    # Clean and process the data
    df = df.dropna(subset=[gene_column, targeting_column])

    if df.empty:
        logger.warning("No valid data found in genomic targeting file after cleanup")
        return {}

    # Convert gene symbols to strings and strip whitespace
    df[gene_column] = df[gene_column].astype(str).str.strip()

    # Convert targeting values to boolean
    targeting_flags = {}

    for _, row in df.iterrows():
        gene_symbol = row[gene_column]
        targeting_value = row[targeting_column]

        # Skip empty gene symbols
        if not gene_symbol or gene_symbol.lower() in ["", "nan", "none"]:
            continue

        # Convert targeting value to boolean
        targeting_bool = convert_to_boolean(targeting_value, default_value)
        targeting_flags[gene_symbol] = targeting_bool

    # Log summary
    true_count = sum(targeting_flags.values())
    false_count = len(targeting_flags) - true_count
    logger.info(
        f"Processed {len(targeting_flags)} genes: {true_count} with targeting=True, "
        f"{false_count} with targeting=False",
    )

    return targeting_flags


def convert_to_boolean(value: Any, default: bool = False) -> bool:
    """
    Convert various value types to boolean for targeting flags.

    Args:
        value: Value to convert (can be string, bool, int, etc.)
        default: Default value if conversion fails

    Returns:
        Boolean representation of the value
    """
    if pd.isna(value):
        return default

    if isinstance(value, bool):
        return value

    if isinstance(value, int | float):
        return bool(value)

    if isinstance(value, str):
        value_lower = value.lower().strip()
        if value_lower in ["true", "yes", "1", "y", "t"]:
            return True
        elif value_lower in ["false", "no", "0", "n", "f"]:
            return False
        else:
            logger.warning(
                f"Unrecognized targeting value '{value}', using default: {default}",
            )
            return default

    logger.warning(
        f"Cannot convert targeting value '{value}' to boolean, using default: {default}",
    )
    return default


def apply_genomic_targeting_flags(
    df: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Apply genomic targeting flags to an annotated gene DataFrame.

    This function adds a 'genomic_targeting' column to the DataFrame based on
    configured targeting flags, with genes not found in the targeting file
    receiving the default value.

    Args:
        df: Annotated DataFrame with gene data
        config: Configuration dictionary

    Returns:
        DataFrame with genomic_targeting column added
    """
    if df.empty:
        logger.info("Empty DataFrame provided for genomic targeting flags")
        return df

    # Check if gene symbol column exists
    if "approved_symbol" not in df.columns:
        logger.error(
            "approved_symbol column not found - cannot apply genomic targeting flags",
        )
        return df

    # Fetch targeting flags
    targeting_flags = fetch_genomic_targeting_flags(config)

    # Get default value from config
    config_manager = ConfigManager(config)
    default_value = config_manager.get_nested(
        "genomic_targeting",
        "default_value",
        default=False,
    )

    # Apply targeting flags
    df_with_targeting = df.copy()
    targeting_values = []

    genes_found = 0
    genes_not_found = 0

    for gene_symbol in df_with_targeting["approved_symbol"]:
        if pd.isna(gene_symbol) or not gene_symbol:
            targeting_values.append(default_value)
            genes_not_found += 1
        elif gene_symbol in targeting_flags:
            targeting_values.append(targeting_flags[gene_symbol])
            genes_found += 1
        else:
            targeting_values.append(default_value)
            genes_not_found += 1

    df_with_targeting["genomic_targeting"] = targeting_values

    # Log summary
    logger.info(
        f"Applied genomic targeting flags: {genes_found} genes found in targeting file, "
        f"{genes_not_found} genes using default value ({default_value})",
    )

    return df_with_targeting


def validate_genomic_targeting_config(config: dict[str, Any]) -> list[str]:
    """
    Validate genomic targeting configuration.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation errors
    """
    errors = []
    config_manager = ConfigManager(config)
    targeting_config = config_manager.get_nested("genomic_targeting", default={})

    if not isinstance(targeting_config, dict):
        errors.append("genomic_targeting config must be a dictionary")
        return errors

    # Check enabled flag
    enabled = targeting_config.get("enabled", False)
    if not isinstance(enabled, bool):
        errors.append("genomic_targeting.enabled must be a boolean")

    # If disabled, skip other validations
    if not enabled:
        return errors

    # Check file_path
    file_path = targeting_config.get("file_path")
    if not file_path:
        errors.append("genomic_targeting.file_path is required when enabled")
    elif not isinstance(file_path, str):
        errors.append("genomic_targeting.file_path must be a string")
    else:
        file_path_obj = Path(file_path)
        allow_missing = targeting_config.get("allow_missing_file", True)
        if not allow_missing and not file_path_obj.exists():
            errors.append(f"genomic_targeting.file_path does not exist: {file_path}")

    # Check column names
    gene_column = targeting_config.get("gene_column", "gene_symbol")
    if not isinstance(gene_column, str):
        errors.append("genomic_targeting.gene_column must be a string")

    targeting_column = targeting_config.get("targeting_column", "targeting")
    if not isinstance(targeting_column, str):
        errors.append("genomic_targeting.targeting_column must be a string")

    # Check default value
    default_value = targeting_config.get("default_value", False)
    if not isinstance(default_value, bool):
        errors.append("genomic_targeting.default_value must be a boolean")

    # Check optional boolean flags
    allow_missing_file = targeting_config.get("allow_missing_file", True)
    if not isinstance(allow_missing_file, bool):
        errors.append("genomic_targeting.allow_missing_file must be a boolean")

    validate_symbols = targeting_config.get("validate_gene_symbols", False)
    if not isinstance(validate_symbols, bool):
        errors.append("genomic_targeting.validate_gene_symbols must be a boolean")

    return errors


def get_genomic_targeting_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary of genomic targeting configuration and data.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    config_manager = ConfigManager(config)
    targeting_config = config_manager.get_nested("genomic_targeting", default={})

    summary = {
        "enabled": targeting_config.get("enabled", False),
        "file_path": targeting_config.get("file_path"),
        "validation_errors": validate_genomic_targeting_config(config),
    }

    # If enabled, try to get flag counts
    if targeting_config.get("enabled", False):
        try:
            targeting_flags = fetch_genomic_targeting_flags(config)
            summary["total_genes"] = len(targeting_flags)
            summary["genes_with_targeting"] = sum(targeting_flags.values())
            summary["genes_without_targeting"] = len(targeting_flags) - sum(
                targeting_flags.values(),
            )
        except Exception as e:
            summary["error"] = str(e)
            summary["total_genes"] = 0
            summary["genes_with_targeting"] = 0
            summary["genes_without_targeting"] = 0
    else:
        summary["total_genes"] = 0
        summary["genes_with_targeting"] = 0
        summary["genes_without_targeting"] = 0

    return summary
