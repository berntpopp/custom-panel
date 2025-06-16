"""
In-house gene panels data source extractor.

This module processes local gene panel files (Excel, CSV, etc.) according to
configuration specifications.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)


def fetch_inhouse_panels_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch gene panel data from in-house panel files.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with in-house panel data
    """
    inhouse_config = config.get("data_sources", {}).get("inhouse_panels", {})

    if not inhouse_config.get("enabled", True):
        logger.info("In-house panels data source is disabled")
        return pd.DataFrame()

    panels_config = inhouse_config.get("panels", [])
    if not panels_config:
        logger.warning("No in-house panels configured")
        return pd.DataFrame()

    all_dataframes = []

    for panel_config in panels_config:
        panel_name = panel_config.get("name", "Unknown_Panel")
        file_path = panel_config.get("file_path")

        if not file_path:
            logger.warning(f"No file path specified for panel {panel_name}")
            continue

        logger.info(f"Processing in-house panel: {panel_name}")

        try:
            df = process_inhouse_panel(panel_config)
            if df is not None and not df.empty:
                all_dataframes.append(df)
                logger.info(f"Successfully processed {len(df)} genes from {panel_name}")
            else:
                logger.warning(f"No genes extracted from {panel_name}")
        except Exception as e:
            logger.error(f"Error processing panel {panel_name}: {e}")
            continue

    # Combine all dataframes
    if all_dataframes:
        combined_df = pd.concat(all_dataframes, ignore_index=True)
        logger.info(
            f"Fetched {len(combined_df)} total gene records from in-house panels"
        )
        return combined_df
    else:
        logger.warning("No in-house panel data was successfully processed")
        return pd.DataFrame()


def process_inhouse_panel(panel_config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Process a single in-house panel file.

    Args:
        panel_config: Panel configuration dictionary

    Returns:
        Standardized DataFrame or None if processing failed
    """
    panel_name = panel_config.get("name", "Unknown_Panel")
    file_path_str = panel_config.get("file_path")
    gene_column = panel_config.get("gene_column", "Gene")
    sheet_name = panel_config.get("sheet_name")
    evidence_score = panel_config.get("evidence_score", 1.0)

    if not file_path_str:
        logger.error(f"No file path specified for panel: {panel_name}")
        return None

    # Validate file path
    file_path = Path(file_path_str)
    if not file_path.exists():
        logger.error(f"File not found: {file_path}")
        return None

    # Read file based on extension
    try:
        if file_path.suffix.lower() in [".xlsx", ".xls"]:
            df = read_excel_panel(file_path, gene_column, sheet_name)
        elif file_path.suffix.lower() == ".csv":
            df = read_csv_panel(file_path, gene_column)
        elif file_path.suffix.lower() == ".tsv":
            df = read_tsv_panel(file_path, gene_column)
        elif file_path.suffix.lower() == ".txt":
            df = read_text_panel(file_path, gene_column)
        else:
            logger.error(f"Unsupported file format: {file_path.suffix}")
            return None
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        return None

    if df is None or df.empty:
        logger.warning(f"No valid data found in {file_path}")
        return None

    # Extract genes from the specified column
    genes = extract_genes_from_column(df, gene_column)

    if not genes:
        logger.warning(f"No genes found in column '{gene_column}' of {file_path}")
        return None

    # Create standardized dataframe
    evidence_scores = [evidence_score] * len(genes)
    source_details = [f"File:{file_path.name}"] * len(genes)

    return create_standard_dataframe(
        genes=genes,
        source_name=panel_name,
        evidence_scores=evidence_scores,
        source_details=source_details,
        gene_names_reported=genes,  # Assume gene names are already standardized
    )


def read_excel_panel(
    file_path: Path, gene_column: str, sheet_name: str | None = None
) -> pd.DataFrame:
    """
    Read gene panel from Excel file.

    Args:
        file_path: Path to Excel file
        gene_column: Name of the gene column
        sheet_name: Name of the sheet to read (optional)

    Returns:
        DataFrame with gene data
    """
    try:
        # Read Excel file
        if sheet_name:
            df = pd.read_excel(file_path, sheet_name=sheet_name, engine="openpyxl")
        else:
            df = pd.read_excel(file_path, engine="openpyxl")

        logger.debug(
            f"Read Excel file {file_path} with {len(df)} rows and columns: {list(df.columns)}"
        )
        return df
    except Exception as e:
        logger.error(f"Error reading Excel file {file_path}: {e}")
        return pd.DataFrame()


def read_csv_panel(file_path: Path, gene_column: str) -> pd.DataFrame:
    """
    Read gene panel from CSV file.

    Args:
        file_path: Path to CSV file
        gene_column: Name of the gene column

    Returns:
        DataFrame with gene data
    """
    try:
        # Try different encodings and separators
        for encoding in ["utf-8", "latin1", "cp1252"]:
            try:
                df = pd.read_csv(file_path, encoding=encoding)
                logger.debug(f"Read CSV file {file_path} with encoding {encoding}")
                return df
            except UnicodeDecodeError:
                continue

        logger.error(f"Could not read CSV file {file_path} with any encoding")
        return pd.DataFrame()
    except Exception as e:
        logger.error(f"Error reading CSV file {file_path}: {e}")
        return pd.DataFrame()


def read_tsv_panel(file_path: Path, gene_column: str) -> pd.DataFrame:
    """
    Read gene panel from TSV file.

    Args:
        file_path: Path to TSV file
        gene_column: Name of the gene column

    Returns:
        DataFrame with gene data
    """
    try:
        # Try different encodings
        for encoding in ["utf-8", "latin1", "cp1252"]:
            try:
                df = pd.read_csv(file_path, sep="\t", encoding=encoding)
                logger.debug(f"Read TSV file {file_path} with encoding {encoding}")
                return df
            except UnicodeDecodeError:
                continue

        logger.error(f"Could not read TSV file {file_path} with any encoding")
        return pd.DataFrame()
    except Exception as e:
        logger.error(f"Error reading TSV file {file_path}: {e}")
        return pd.DataFrame()


def read_text_panel(file_path: Path, gene_column: str) -> pd.DataFrame:
    """
    Read gene panel from plain text file (one gene per line).

    Args:
        file_path: Path to text file
        gene_column: Name of the gene column (used as column name)

    Returns:
        DataFrame with gene data
    """
    try:
        with open(file_path, encoding="utf-8") as f:
            genes = [line.strip() for line in f if line.strip()]

        if not genes:
            return pd.DataFrame()

        # Create DataFrame with genes
        df = pd.DataFrame({gene_column: genes})
        logger.debug(f"Read text file {file_path} with {len(genes)} genes")
        return df
    except Exception as e:
        logger.error(f"Error reading text file {file_path}: {e}")
        return pd.DataFrame()


def extract_genes_from_column(df: pd.DataFrame, gene_column: str) -> list[str]:
    """
    Extract gene symbols from a DataFrame column.

    Args:
        df: Input DataFrame
        gene_column: Name of the column containing gene symbols

    Returns:
        List of unique gene symbols
    """
    if gene_column not in df.columns:
        logger.error(
            f"Column '{gene_column}' not found in DataFrame. Available columns: {list(df.columns)}"
        )
        return []

    # Extract genes from the column
    genes_series = df[gene_column].dropna()

    # Handle different data types
    genes = []
    for gene in genes_series:
        if isinstance(gene, str):
            # Handle comma-separated genes in a single cell
            if "," in gene:
                genes.extend([g.strip() for g in gene.split(",") if g.strip()])
            elif ";" in gene:
                genes.extend([g.strip() for g in gene.split(";") if g.strip()])
            else:
                genes.append(gene.strip())
        elif pd.notna(gene):
            # Convert non-string values to string
            genes.append(str(gene).strip())

    # Remove empty strings and duplicates while preserving order
    unique_genes = []
    seen = set()
    for gene in genes:
        if gene and gene not in seen:
            unique_genes.append(gene)
            seen.add(gene)

    logger.debug(
        f"Extracted {len(unique_genes)} unique genes from column '{gene_column}'"
    )
    return unique_genes


def validate_inhouse_panel_config(panel_config: dict[str, Any]) -> list[str]:
    """
    Validate in-house panel configuration.

    Args:
        panel_config: Panel configuration dictionary

    Returns:
        List of validation errors (empty if valid)
    """
    errors = []

    # Check required fields
    if not panel_config.get("name"):
        errors.append("Panel name is required")

    if not panel_config.get("file_path"):
        errors.append("File path is required")

    if not panel_config.get("gene_column"):
        errors.append("Gene column name is required")

    # Check file existence
    file_path = panel_config.get("file_path")
    if file_path and not Path(file_path).exists():
        errors.append(f"File not found: {file_path}")

    # Check evidence score
    evidence_score = panel_config.get("evidence_score", 1.0)
    if not isinstance(evidence_score, int | float) or evidence_score < 0:
        errors.append("Evidence score must be a non-negative number")

    return errors


def get_inhouse_panel_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary information about configured in-house panels.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    inhouse_config = config.get("data_sources", {}).get("inhouse_panels", {})
    panels_config = inhouse_config.get("panels", [])

    summary = {
        "enabled": inhouse_config.get("enabled", True),
        "total_panels": len(panels_config),
        "panels": [],
    }

    for panel_config in panels_config:
        panel_summary = {
            "name": panel_config.get("name", "Unknown"),
            "file_path": panel_config.get("file_path", ""),
            "gene_column": panel_config.get("gene_column", ""),
            "evidence_score": panel_config.get("evidence_score", 1.0),
            "file_exists": Path(panel_config.get("file_path", "")).exists()
            if panel_config.get("file_path")
            else False,
            "validation_errors": validate_inhouse_panel_config(panel_config),
        }
        summary["panels"].append(panel_summary)

    return summary
