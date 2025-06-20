"""
Polygenic Risk Score (PRS) SNPs source fetcher.

This module fetches PRS SNPs from configured panel files.
These SNPs are used in risk score calculations and typically
contain additional metadata like effect alleles and weights.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..parsers_snp.parsers_prs import create_prs_parser

logger = logging.getLogger(__name__)


def fetch_prs_snps(config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch PRS SNPs from configured panels.

    Args:
        config: Configuration dictionary

    Returns:
        DataFrame with PRS SNPs or None if disabled/failed
    """
    snp_config = config.get("snp_processing", {})
    
    if not snp_config.get("enabled", False):
        logger.info("SNP processing is disabled")
        return None

    prs_config = snp_config.get("prs", {})
    
    if not prs_config.get("enabled", False):
        logger.info("PRS SNPs are disabled")
        return None

    panels = prs_config.get("panels", [])
    if not panels:
        logger.warning("No PRS panels configured")
        return None

    logger.info(f"Fetching PRS SNPs from {len(panels)} panels")

    all_snps = []

    for panel_config in panels:
        try:
            panel_df = _fetch_single_prs_panel(panel_config)
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
            logger.error(
                f"✗ {panel_config.get('name', 'Unknown')}: {e}"
            )

    if not all_snps:
        logger.warning("No PRS SNPs were successfully fetched")
        return None

    # Combine all panels
    combined_df = pd.concat(all_snps, ignore_index=True)
    
    # Check for conflicting effect alleles before aggregation
    if "effect_allele" in combined_df.columns:
        conflicts = _check_effect_allele_conflicts(combined_df)
        if conflicts:
            logger.warning(f"Found {len(conflicts)} rsIDs with conflicting effect alleles")
            for rsid, alleles in conflicts.items():
                logger.warning(f"  {rsid}: {alleles}")
    
    # Apply R-script-like aggregation but preserve PRS metadata
    prs_snps_panel = _aggregate_prs_snps_by_rsid(combined_df)

    logger.info(f"Successfully fetched {len(prs_snps_panel)} unique PRS SNPs")
    return prs_snps_panel


def _fetch_single_prs_panel(panel_config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch SNPs from a single PRS panel.

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
        parser = create_prs_parser(file_path, panel_config)
        df = parser.parse()
        
        if df.empty:
            logger.warning(f"Panel {name} contains no valid SNPs")
            return None
            
        return df
        
    except Exception as e:
        raise Exception(f"Failed to parse panel {name}: {e}") from e


def _aggregate_prs_snps_by_rsid(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate PRS SNPs by rsID, merging sources while preserving metadata.
    
    This function replicates the R script aggregation pattern but preserves
    PRS-specific metadata like effect alleles, weights, etc.

    Args:
        df: DataFrame with individual PRS SNP records

    Returns:
        DataFrame with aggregated PRS SNPs
    """
    # Rename rsid to snp to match R script column names
    if "rsid" in df.columns:
        df = df.rename(columns={"rsid": "snp"})
    
    # Define aggregation rules for different column types
    agg_rules = {
        "source": lambda x: "; ".join(x.dropna().astype(str).unique()),
        "category": "first"
    }
    
    # For PRS metadata, keep first non-null value
    prs_metadata_cols = ["chromosome", "position", "effect_allele", "effect_weight", 
                        "OR", "beta", "se", "pvalue", "freq"]
    
    for col in prs_metadata_cols:
        if col in df.columns:
            agg_rules[col] = "first"  # Keep first non-null value
    
    # Group by snp and aggregate
    aggregated = df.groupby("snp").agg(agg_rules).reset_index()
    
    # Sort by snp (matching R script: arrange(snp))
    aggregated = aggregated.sort_values("snp").reset_index(drop=True)
    
    return aggregated


def _check_effect_allele_conflicts(df: pd.DataFrame) -> dict[str, list[str]]:
    """
    Check for rsIDs with conflicting effect alleles across panels.

    Args:
        df: DataFrame with PRS SNPs

    Returns:
        Dictionary mapping rsIDs to list of conflicting effect alleles
    """
    if "effect_allele" not in df.columns:
        return {}

    conflicts = {}
    
    # Use snp column if rsid was renamed
    rsid_col = "snp" if "snp" in df.columns else "rsid"
    
    # Group by rsID and check for multiple effect alleles
    for rsid, group in df.groupby(rsid_col):
        unique_alleles = group["effect_allele"].dropna().unique()
        if len(unique_alleles) > 1:
            conflicts[rsid] = unique_alleles.tolist()
    
    return conflicts


def get_prs_snps_summary(df: pd.DataFrame) -> dict[str, Any]:
    """
    Generate summary statistics for PRS SNPs.

    Args:
        df: DataFrame with PRS SNPs

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

    # Effect allele information
    if "effect_allele" in df.columns:
        summary["with_effect_allele"] = df["effect_allele"].notna().sum()
        summary["effect_allele_distribution"] = df["effect_allele"].value_counts().to_dict()

    # Chromosome distribution
    if "chromosome" in df.columns:
        summary["chromosome_distribution"] = df["chromosome"].value_counts().to_dict()

    # Additional PRS-specific metrics
    prs_columns = ["effect_weight", "OR", "beta", "pvalue"]
    for col in prs_columns:
        if col in df.columns:
            summary[f"with_{col}"] = df[col].notna().sum()

    return summary