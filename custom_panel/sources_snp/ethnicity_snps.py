"""
Ethnicity SNPs source fetcher.

This module fetches ethnicity SNPs from configured files, specifically
handling the ethnicity_snps.xlsx file format used in the original R implementation.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)


def fetch_ethnicity_snps(config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch ethnicity SNPs from configured ethnicity files.

    This function specifically handles the ethnicity_snps.xlsx format
    from the original R implementation (V03_Ethnicity). Returns raw,
    unharmonized SNP data. Harmonization is handled centrally in the
    Pipeline class for improved efficiency.

    Args:
        config: Configuration dictionary

    Returns:
        DataFrame with raw ethnicity SNPs or None if disabled/failed
    """
    snp_config = config.get("snp_processing", {})

    if not snp_config.get("enabled", False):
        logger.info("SNP processing is disabled")
        return None

    ethnicity_config = snp_config.get("ethnicity", {})

    if not ethnicity_config.get("enabled", False):
        logger.info("Ethnicity SNPs are disabled")
        return None

    # Get ethnicity panels from configuration
    panels = ethnicity_config.get("panels", [])

    if not panels:
        logger.warning("No ethnicity panels configured")
        return None

    logger.info(f"Fetching ethnicity SNPs from {len(panels)} panels")

    all_snps = []

    for panel_config in panels:
        try:
            panel_df = _fetch_single_ethnicity_panel(panel_config)
            if panel_df is not None and not panel_df.empty:
                all_snps.append(panel_df)
                logger.info(
                    f"✓ {panel_config.get('name', 'Unknown')}: {len(panel_df)} SNPs"
                )
            else:
                logger.warning(
                    f"⚠ {panel_config.get('name', 'Unknown')}: No SNPs found"
                )
        except Exception as e:
            logger.error(f"✗ {panel_config.get('name', 'Unknown')}: {e}")

    if not all_snps:
        logger.warning("No ethnicity SNPs were successfully fetched")
        return None

    # Combine all panels and apply R-script-like aggregation
    combined_df = pd.concat(all_snps, ignore_index=True)

    # Group by rsID and merge sources with "; " separator (matching R implementation)
    ethnicity_snps_panel = _aggregate_snps_by_rsid(combined_df)

    logger.info(
        f"Successfully fetched {len(ethnicity_snps_panel)} unique ethnicity SNPs"
    )
    return ethnicity_snps_panel


def _fetch_single_ethnicity_panel(panel_config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch SNPs from a single ethnicity panel (Excel format).

    Args:
        panel_config: Configuration for a single panel

    Returns:
        DataFrame with SNPs from the panel or None if failed

    Raises:
        Exception: If panel fetching fails
    """
    name = panel_config.get("name", "Unknown")
    file_path = panel_config.get("file_path")

    if not file_path:
        raise ValueError(f"No file_path specified for panel {name}")

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Panel file not found: {file_path}")

    # Get column mappings from config
    rsid_column = panel_config.get("rsid_column", "rs_id")
    source_column = panel_config.get("source_column", "group")
    sheet_name = panel_config.get("sheet_name", 0)  # Default to first sheet

    try:
        # Read Excel file (similar to read_excel in R)
        df = pd.read_excel(file_path, sheet_name=sheet_name)

        # Check for required columns
        if rsid_column not in df.columns:
            raise ValueError(
                f"Required rsID column '{rsid_column}' not found. Available columns: {list(df.columns)}"
            )

        # Extract and validate data
        if source_column in df.columns:
            # We have both rsID and source group columns
            result_df = df[[rsid_column, source_column]].copy()
            result_df.columns = ["snp", "source_group"]

            # Remove rows with missing rsIDs
            result_df = result_df.dropna(subset=["snp"])

            # Add standard columns
            result_df["source"] = name
            result_df["category"] = "ethnicity"

        else:
            # Only rsID column available
            result_df = pd.DataFrame(
                {
                    "snp": df[rsid_column].dropna(),
                    "source": name,
                    "category": "ethnicity",
                }
            )

        if result_df.empty:
            logger.warning(f"Panel {name} contains no valid SNPs")
            return None

        # Standardize rsID format
        result_df["snp"] = result_df["snp"].apply(_standardize_rsid)

        # Remove any rows with empty rsids
        result_df = result_df.dropna(subset=["snp"])
        result_df = result_df[result_df["snp"].str.strip() != ""]

        return result_df

    except Exception as e:
        raise Exception(f"Failed to parse ethnicity panel {name}: {e}") from e


def _aggregate_snps_by_rsid(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate SNPs by rsID, merging sources with "; " separator.

    This function replicates the R script pattern:
    ```r
    group_by(rs_id) %>%
    summarise(source = paste0(group, collapse = "; "))
    ```

    Args:
        df: DataFrame with individual SNP records

    Returns:
        DataFrame with aggregated SNPs
    """
    if "source_group" in df.columns:
        # Use source_group for aggregation (matching R script)
        aggregated = (
            df.groupby("snp")
            .agg(
                {
                    "source_group": lambda x: "; ".join(
                        x.dropna().astype(str).unique()
                    ),
                    "source": "first",  # Keep first source name
                    "category": "first",  # Keep first category
                }
            )
            .reset_index()
        )

        # Rename source_group to source (matching R output)
        aggregated = aggregated.rename(columns={"source_group": "source_detail"})
        aggregated["source"] = aggregated[
            "source_detail"
        ]  # Use aggregated source_groups as main source

    else:
        # Standard aggregation by source name
        aggregated = (
            df.groupby("snp")
            .agg(
                {
                    "source": lambda x: "; ".join(x.dropna().astype(str).unique()),
                    "category": "first",
                }
            )
            .reset_index()
        )

    # Sort by snp (matching R script: arrange(snp))
    aggregated = aggregated.sort_values("snp").reset_index(drop=True)

    return aggregated


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


def get_ethnicity_snps_summary(df: pd.DataFrame) -> dict[str, Any]:
    """
    Generate summary statistics for ethnicity SNPs.

    Args:
        df: DataFrame with ethnicity SNPs

    Returns:
        Summary statistics dictionary
    """
    if df.empty:
        return {"total_snps": 0}

    summary = {
        "total_snps": len(df),
        "unique_rsids": df["snp"].nunique(),
        "categories": df["category"].value_counts().to_dict(),
    }

    # Source breakdown
    if "source" in df.columns:
        summary["source_breakdown"] = df["source"].value_counts().to_dict()

    # Source detail breakdown if available
    if "source_detail" in df.columns:
        summary["source_detail_breakdown"] = (
            df["source_detail"].value_counts().to_dict()
        )

    return summary
