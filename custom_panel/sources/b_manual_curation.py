"""
Manual curation data source extractor.

This module processes manually curated gene lists from various file formats
(Excel, CSV, TXT) according to configuration specifications.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.io import create_standard_dataframe
from .g00_inhouse_panels import (
    extract_genes_from_column,
    read_csv_panel,
    read_excel_panel,
    read_text_panel,
)

logger = logging.getLogger(__name__)


def fetch_manual_curation_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch manually curated gene data from configured file lists.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with manual curation data
    """
    manual_config = config.get("data_sources", {}).get("Manual_Curation", {})

    if not manual_config.get("enabled", True):
        logger.info("Manual curation data source is disabled")
        return pd.DataFrame()

    lists_config = manual_config.get("lists", [])
    if not lists_config:
        logger.warning("No manual curation lists configured")
        return pd.DataFrame()

    all_dataframes = []

    for list_config in lists_config:
        try:
            df = process_manual_list(list_config)
            if not df.empty:
                all_dataframes.append(df)
        except Exception as e:
            logger.error(
                f"Failed to process manual list {list_config.get('name', 'unknown')}: {e}"
            )
            continue

    if not all_dataframes:
        logger.warning("No manual curation data successfully processed")
        return pd.DataFrame()

    # Combine all dataframes
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    logger.info(
        f"Successfully processed {len(combined_df)} genes from {len(all_dataframes)} manual lists"
    )

    return combined_df


def process_manual_list(list_config: dict[str, Any]) -> pd.DataFrame:
    """
    Process a single manual curation list.

    Args:
        list_config: Configuration for a single manual list

    Returns:
        Standardized DataFrame for the list
    """
    name = list_config.get("name", "Unknown")
    file_path = list_config.get("file_path")
    gene_column = list_config.get("gene_column", "gene_symbol")
    evidence_score = list_config.get("evidence_score", 1.0)
    category = list_config.get("category", "germline")  # Default to germline

    if not file_path:
        logger.error(f"No file_path specified for manual list '{name}'")
        return pd.DataFrame()

    file_path = Path(file_path)
    if not file_path.exists():
        logger.error(f"File not found for manual list '{name}': {file_path}")
        return pd.DataFrame()

    logger.info(f"Processing manual list '{name}' from {file_path}")

    # Read the file based on extension
    try:
        file_extension = file_path.suffix.lower()

        if file_extension in [".xlsx", ".xls"]:
            df = read_excel_panel(file_path, gene_column)
        elif file_extension == ".csv":
            df = read_csv_panel(file_path, gene_column)
        elif file_extension in [".txt", ".tsv"]:
            df = read_text_panel(file_path, gene_column)
        else:
            logger.error(
                f"Unsupported file format '{file_extension}' for manual list '{name}'"
            )
            return pd.DataFrame()

        if df.empty:
            logger.warning(f"No data found in file for manual list '{name}'")
            return pd.DataFrame()

    except Exception as e:
        logger.error(f"Error reading file for manual list '{name}': {e}")
        return pd.DataFrame()

    # Extract genes from the specified column
    try:
        genes = extract_genes_from_column(df, gene_column)
        if not genes:
            logger.warning(
                f"No genes extracted from column '{gene_column}' in manual list '{name}'"
            )
            return pd.DataFrame()

        logger.info(f"Extracted {len(genes)} genes from manual list '{name}'")

    except Exception as e:
        logger.error(f"Error extracting genes from manual list '{name}': {e}")
        return pd.DataFrame()

    # Create standardized dataframe
    source_name = f"Manual_Curation:{name}"
    evidence_scores = [evidence_score] * len(genes)
    source_details = [
        f"File:{file_path.name}|Column:{gene_column}|Category:{category}"
    ] * len(genes)

    standardized_df = create_standard_dataframe(
        genes=genes,
        source_name=source_name,
        evidence_scores=evidence_scores,
        source_details=source_details,
        gene_names_reported=genes,
    )

    # Add temporary category column for use by the merger
    standardized_df["category"] = category

    logger.info(
        f"Created standardized dataframe for manual list '{name}' with {len(standardized_df)} genes (category: {category})"
    )
    return standardized_df


def validate_manual_curation_config(config: dict[str, Any]) -> list[str]:
    """
    Validate manual curation configuration.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation errors
    """
    errors = []
    manual_config = config.get("data_sources", {}).get("Manual_Curation", {})

    if not isinstance(manual_config, dict):
        errors.append("manual_curation config must be a dictionary")
        return errors

    lists_config = manual_config.get("lists", [])
    if not isinstance(lists_config, list):
        errors.append("manual_curation.lists must be a list")
        return errors

    for i, list_config in enumerate(lists_config):
        if not isinstance(list_config, dict):
            errors.append(f"manual_curation.lists[{i}] must be a dictionary")
            continue

        # Check required fields
        name = list_config.get("name")
        if not name or not isinstance(name, str):
            errors.append(
                f"manual_curation.lists[{i}].name is required and must be a string"
            )

        file_path = list_config.get("file_path")
        if not file_path or not isinstance(file_path, str):
            errors.append(
                f"manual_curation.lists[{i}].file_path is required and must be a string"
            )
        elif not Path(file_path).exists():
            errors.append(
                f"manual_curation.lists[{i}].file_path does not exist: {file_path}"
            )

        # Check optional fields
        evidence_score = list_config.get("evidence_score", 1.0)
        if (
            not isinstance(evidence_score, int | float)
            or evidence_score < 0
            or evidence_score > 1
        ):
            errors.append(
                f"manual_curation.lists[{i}].evidence_score must be a number between 0 and 1"
            )

        gene_column = list_config.get("gene_column", "gene_symbol")
        if not isinstance(gene_column, str):
            errors.append(f"manual_curation.lists[{i}].gene_column must be a string")

    return errors


def get_manual_curation_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary of manual curation configuration.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    manual_config = config.get("data_sources", {}).get("Manual_Curation", {})

    summary = {
        "enabled": manual_config.get("enabled", True),
        "lists_count": len(manual_config.get("lists", [])),
        "validation_errors": validate_manual_curation_config(config),
    }

    # Count total genes across all lists
    total_genes = 0
    lists_config = manual_config.get("lists", [])
    for list_config in lists_config:
        try:
            df = process_manual_list(list_config)
            total_genes += len(df)
        except Exception:
            continue

    summary["total_genes"] = total_genes
    return summary
