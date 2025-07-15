"""
Manual SNPs source fetcher.

This module fetches manually curated SNPs from configured lists.
These SNPs are typically added for specific purposes like pharmacogenomics (PGx).
"""

import logging
from pathlib import Path
from typing import Any, Optional

import pandas as pd

logger = logging.getLogger(__name__)


def fetch_manual_snps(
    config: dict[str, Any], harmonizer: Optional[Any] = None
) -> pd.DataFrame | None:
    """
    Fetch manually curated SNPs from configured lists.

    Args:
        config: Configuration dictionary

    Returns:
        DataFrame with manual SNPs or None if disabled/failed
    """
    snp_config = config.get("snp_processing", {})

    if not snp_config.get("enabled", False):
        logger.info("SNP processing is disabled")
        return None

    manual_config = snp_config.get("manual_snps", {})

    if not manual_config.get("enabled", False):
        logger.info("Manual SNPs are disabled")
        return None

    lists = manual_config.get("lists", [])
    if not lists:
        logger.warning("No manual SNP lists configured")
        return None

    logger.info(f"Fetching manual SNPs from {len(lists)} lists")

    all_snps = []

    for list_config in lists:
        try:
            list_df = _fetch_single_manual_list(list_config)
            if list_df is not None and not list_df.empty:
                all_snps.append(list_df)
                logger.info(
                    f"✓ {list_config.get('name', 'Unknown')}: {len(list_df)} SNPs"
                )
            else:
                logger.warning(f"⚠ {list_config.get('name', 'Unknown')}: No SNPs found")
        except Exception as e:
            logger.error(f"✗ {list_config.get('name', 'Unknown')}: {e}")

    if not all_snps:
        logger.warning("No manual SNPs were successfully fetched")
        return None

    # Combine all lists
    combined_df = pd.concat(all_snps, ignore_index=True)

    # Remove duplicates based on SNP ID, keeping first occurrence
    initial_count = len(combined_df)
    combined_df = combined_df.drop_duplicates(subset=["snp"], keep="first")
    final_count = len(combined_df)

    if initial_count != final_count:
        logger.info(
            f"Removed {initial_count - final_count} duplicate rsIDs from manual lists"
        )

    # Apply harmonization if harmonizer is provided
    if harmonizer is not None:
        try:
            logger.info(f"Harmonizing {len(combined_df)} manual SNPs")
            harmonized_manual_snps = harmonizer.harmonize_snp_batch(combined_df)

            if not harmonized_manual_snps.empty:
                logger.info(
                    f"Successfully harmonized {len(harmonized_manual_snps)} manual SNPs"
                )
                return harmonized_manual_snps
            else:
                logger.warning(
                    "Harmonization resulted in empty DataFrame, returning original data"
                )

        except Exception as e:
            logger.error(f"Error during manual SNP harmonization: {e}")
            logger.info("Continuing with non-harmonized data")

    logger.info(f"Successfully fetched {final_count} unique manual SNPs")
    return combined_df


def _fetch_single_manual_list(list_config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch SNPs from a single manual list.

    Args:
        list_config: Configuration for a single manual list

    Returns:
        DataFrame with SNPs from the list or None if failed

    Raises:
        Exception: If list fetching fails
    """
    name = list_config.get("name", "Unknown")
    file_path = list_config.get("file_path")

    if not file_path:
        raise ValueError(f"No file_path specified for manual list {name}")

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Manual list file not found: {file_path}")

    # Get parser configuration
    parser_type = list_config.get("parser", "manual_csv")
    rsid_column = list_config.get("rsid_column", "rsID")
    details_column = list_config.get("details_column")

    try:
        if parser_type == "manual_csv":
            df = _parse_manual_csv(file_path, rsid_column, details_column, name)
        elif parser_type == "manual_tsv":
            df = _parse_manual_tsv(file_path, rsid_column, details_column, name)
        elif parser_type == "manual_excel":
            sheet_name = list_config.get("sheet_name", 0)
            df = _parse_manual_excel(
                file_path, rsid_column, details_column, name, sheet_name
            )
        else:
            raise ValueError(f"Unsupported manual parser type: {parser_type}")

        if df.empty:
            logger.warning(f"Manual list {name} contains no valid SNPs")
            return None

        return df

    except Exception as e:
        raise Exception(f"Failed to parse manual list {name}: {e}") from e


def _parse_manual_csv(
    file_path: Path, rsid_column: str, details_column: str | None, name: str
) -> pd.DataFrame:
    """Parse manual CSV file."""
    df = pd.read_csv(file_path)
    return _process_manual_dataframe(df, rsid_column, details_column, name)


def _parse_manual_tsv(
    file_path: Path, rsid_column: str, details_column: str | None, name: str
) -> pd.DataFrame:
    """Parse manual TSV file."""
    df = pd.read_csv(file_path, sep="\t")
    return _process_manual_dataframe(df, rsid_column, details_column, name)


def _parse_manual_excel(
    file_path: Path,
    rsid_column: str,
    details_column: str | None,
    name: str,
    sheet_name: str | int = 0,
) -> pd.DataFrame:
    """Parse manual Excel file."""
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    return _process_manual_dataframe(df, rsid_column, details_column, name)


def _process_manual_dataframe(
    df: pd.DataFrame, rsid_column: str, details_column: str | None, name: str
) -> pd.DataFrame:
    """
    Process a manual DataFrame to extract SNPs.

    Args:
        df: Input DataFrame
        rsid_column: Column containing rsIDs
        details_column: Optional column with additional details
        name: Name of the manual list

    Returns:
        Processed DataFrame with standardized format
    """
    # Check for required rsID column
    if rsid_column not in df.columns:
        raise ValueError(
            f"Required rsID column '{rsid_column}' not found. "
            f"Available columns: {list(df.columns)}"
        )

    # Extract rsIDs
    rsids = df[rsid_column].dropna().tolist()

    if not rsids:
        raise ValueError(f"No valid rsIDs found in column '{rsid_column}'")

    # Create result DataFrame with 'snp' column for harmonizer compatibility
    result_df = pd.DataFrame(
        {
            "snp": rsids,
            "rsid": rsids,
            "source": name,
            "category": "manual",
            "snp_type": "manual",
        }
    )

    # Add details column if specified and available
    if details_column and details_column in df.columns:
        details = df[details_column].dropna().tolist()
        if len(details) == len(rsids):
            result_df["details"] = details
        else:
            logger.warning(
                f"Details column '{details_column}' length doesn't match rsID column"
            )

    # Standardize rsID format for both columns
    result_df["snp"] = result_df["snp"].apply(_standardize_rsid)
    result_df["rsid"] = result_df["rsid"].apply(_standardize_rsid)

    # Remove any rows with empty rsids
    result_df = result_df.dropna(subset=["snp"])
    result_df = result_df[result_df["snp"].str.strip() != ""]

    return result_df


def _standardize_rsid(rsid: str | None) -> str | None:
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

    # For other formats, return as-is
    return rsid


def get_manual_snps_summary(df: pd.DataFrame) -> dict[str, Any]:
    """
    Generate summary statistics for manual SNPs.

    Args:
        df: DataFrame with manual SNPs

    Returns:
        Summary statistics dictionary
    """
    if df.empty:
        return {"total_snps": 0}

    summary = {
        "total_snps": len(df),
        "unique_rsids": df["rsid"].nunique(),
        "sources": df["source"].nunique(),
        "categories": df["category"].value_counts().to_dict(),
    }

    # Source breakdown
    if "source" in df.columns:
        summary["source_breakdown"] = df["source"].value_counts().to_dict()

    # Details information
    if "details" in df.columns:
        summary["with_details"] = df["details"].notna().sum()
        summary["details_breakdown"] = df["details"].value_counts().to_dict()

    return summary
