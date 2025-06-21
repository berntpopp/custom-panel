"""
Identity SNPs source fetcher.

This module fetches identity SNPs from configured panel files.
These SNPs are typically used for sample tracking and identification.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..parsers_snp.parsers_identity import create_identity_parser

logger = logging.getLogger(__name__)


def fetch_identity_snps(config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch identity SNPs from configured panels.

    Args:
        config: Configuration dictionary

    Returns:
        DataFrame with identity SNPs or None if disabled/failed
    """
    snp_config = config.get("snp_processing", {})

    if not snp_config.get("enabled", False):
        logger.info("SNP processing is disabled")
        return None

    identity_config = snp_config.get("identity", {})

    if not identity_config.get("enabled", False):
        logger.info("Identity SNPs are disabled")
        return None

    panels = identity_config.get("panels", [])
    if not panels:
        logger.warning("No identity panels configured")
        return None

    logger.info(f"Fetching identity SNPs from {len(panels)} panels")

    all_snps = []

    for panel_config in panels:
        try:
            panel_df = _fetch_single_identity_panel(panel_config)
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
        logger.warning("No identity SNPs were successfully fetched")
        return None

    # Combine all panels and apply R-script-like aggregation
    combined_df = pd.concat(all_snps, ignore_index=True)

    # Group by rsID and merge sources with "; " separator (matching R implementation)
    identity_snps_panel = _aggregate_snps_by_rsid(combined_df)

    logger.info(
        f"Successfully fetched {len(identity_snps_panel)} unique identity SNPs"
    )
    return identity_snps_panel


def _fetch_single_identity_panel(panel_config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch SNPs from a single identity panel.

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

    # Create the appropriate parser
    try:
        parser = create_identity_parser(file_path, panel_config)
        df = parser.parse()

        if df.empty:
            logger.warning(f"Panel {name} contains no valid SNPs")
            return None

        return df

    except Exception as e:
        raise Exception(f"Failed to parse panel {name}: {e}") from e


def _aggregate_snps_by_rsid(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate SNPs by rsID, merging sources with "; " separator.

    This function replicates the R script pattern:
    ```r
    group_by(snp) %>%
    summarise(source = paste0(panel_name, collapse = "; "))
    ```

    Args:
        df: DataFrame with individual SNP records

    Returns:
        DataFrame with aggregated SNPs
    """
    # Rename rsid to snp to match R script column names
    if "rsid" in df.columns:
        df = df.rename(columns={"rsid": "snp"})

    # Group by snp and aggregate sources (matching R script)
    aggregated = (
        df.groupby("snp")
        .agg(
            {
                "source": lambda x: "; ".join(x.dropna().astype(str).unique()),
                "category": "first",  # Keep first category
            }
        )
        .reset_index()
    )

    # Sort by snp (matching R script: arrange(snp))
    aggregated = aggregated.sort_values("snp").reset_index(drop=True)

    return aggregated


def get_identity_snps_summary(df: pd.DataFrame) -> dict[str, Any]:
    """
    Generate summary statistics for identity/ethnicity SNPs.

    Args:
        df: DataFrame with identity/ethnicity SNPs

    Returns:
        Summary statistics dictionary
    """
    if df.empty:
        return {"total_snps": 0}

    # Handle both rsid and snp column names
    snp_col = "snp" if "snp" in df.columns else "rsid"

    summary = {
        "total_snps": len(df),
        "unique_rsids": df[snp_col].nunique(),
        "categories": df["category"].value_counts().to_dict(),
    }

    # Source breakdown
    if "source" in df.columns:
        summary["source_breakdown"] = df["source"].value_counts().to_dict()

    # Source group breakdown if available
    if "source_group" in df.columns:
        summary["source_group_breakdown"] = df["source_group"].value_counts().to_dict()

    return summary
