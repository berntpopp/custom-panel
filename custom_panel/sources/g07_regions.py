"""
Regions data source extractor.

This module processes genomic regions from multiple sources:
1. Manual regions from Excel files (data/regions/regions_list.xlsx)
2. Stripy repeat regions from JSON files (data/regions/stripy_variant_catalog.json)
"""

import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.config_manager import ConfigManager

logger = logging.getLogger(__name__)


def fetch_regions_data(config: dict[str, Any]) -> dict[str, pd.DataFrame]:
    """
    Fetch genomic regions data from configured sources.

    Args:
        config: Configuration dictionary

    Returns:
        Dictionary of DataFrames by region source type
    """
    config_manager = ConfigManager(config)
    regions_config = config_manager.to_dict().get("regions_processing", {})

    if not regions_config.get("enabled", False):
        logger.info("Regions processing is disabled")
        return {}

    regions_data = {}

    # Load manual regions if enabled
    manual_config = regions_config.get("manual_regions", {})
    if manual_config.get("enabled", False):
        logger.info("Loading manual regions from Excel file")
        manual_df = _load_manual_regions(manual_config)
        if not manual_df.empty:
            regions_data["manual"] = manual_df
            logger.info(f"Loaded {len(manual_df)} manual regions")

    # Load stripy regions if enabled
    stripy_config = regions_config.get("stripy_regions", {})
    if stripy_config.get("enabled", False):
        logger.info("Loading stripy repeat regions from JSON file")
        stripy_df = _load_stripy_regions(stripy_config)
        if not stripy_df.empty:
            regions_data["stripy"] = stripy_df
            logger.info(f"Loaded {len(stripy_df)} stripy repeat regions")

    logger.info(f"Total regions sources loaded: {len(regions_data)}")
    return regions_data


def _load_manual_regions(config: dict[str, Any]) -> pd.DataFrame:
    """
    Load manual regions from Excel file.

    Args:
        config: Manual regions configuration

    Returns:
        Standardized DataFrame with manual regions
    """
    file_path = Path(config.get("file_path", ""))
    if not file_path.exists():
        logger.error(f"Manual regions file not found: {file_path}")
        return pd.DataFrame()

    try:
        # Read Excel file
        sheet_name = config.get("sheet_name", 0)  # Default to first sheet
        df = pd.read_excel(file_path, sheet_name=sheet_name)

        logger.info(f"Read {len(df)} rows from manual regions file")

        # Filter for included regions
        filter_column = config.get("filter_column", "include")
        filter_value = config.get("filter_value", "yes")

        if filter_column in df.columns:
            df = df[df[filter_column].astype(str).str.lower() == filter_value.lower()]
            logger.info(
                f"Filtered to {len(df)} regions with {filter_column}={filter_value}"
            )

        # Clean column names (remove trailing spaces)
        df.columns = df.columns.str.strip()

        # Check required columns
        required_cols = ["hg38_chr", "hg38_start", "hg38_stop"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            return pd.DataFrame()

        # Standardize the data
        standardized_data: list[dict[str, Any]] = []
        for _, row in df.iterrows():
            region_name = f"manual_{row.get('number', len(standardized_data) + 1)}"

            # Clean chromosome name
            chromosome = str(row["hg38_chr"]).replace("chr", "")

            # Convert coordinates to integers
            try:
                start = int(row["hg38_start"])
                end = int(row["hg38_stop"])
            except (ValueError, TypeError):
                logger.warning(
                    f"Invalid coordinates for region {region_name}, skipping"
                )
                continue

            comment = str(row.get("comment", "")).strip()

            # Calculate coverage in base pairs
            coverage = end - start + 1

            # Get type column, default to "complete" if not specified
            region_type = str(row.get("type", "complete")).strip()

            standardized_data.append(
                {
                    "region_name": region_name,
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "coverage": coverage,
                    "source_type": "manual",
                    "region_type": region_type,
                    "comment": comment,
                }
            )

        if not standardized_data:
            logger.warning("No valid manual regions found")
            return pd.DataFrame()

        result_df = pd.DataFrame(standardized_data)
        logger.info(f"Successfully processed {len(result_df)} manual regions")
        return result_df

    except Exception as e:
        logger.error(f"Error loading manual regions: {e}")
        return pd.DataFrame()


def _load_stripy_regions(config: dict[str, Any]) -> pd.DataFrame:
    """
    Load stripy repeat regions from JSON file.

    Args:
        config: Stripy regions configuration

    Returns:
        Standardized DataFrame with stripy repeat regions
    """
    file_path = Path(config.get("file_path", ""))
    if not file_path.exists():
        logger.error(f"Stripy regions file not found: {file_path}")
        return pd.DataFrame()

    try:
        # Read JSON file
        with open(file_path, encoding="utf-8") as f:
            data = json.load(f)

        if not isinstance(data, list):
            logger.error("Expected JSON array format for stripy regions")
            return pd.DataFrame()

        logger.info(f"Read {len(data)} stripy repeat regions from JSON")

        # Standardize the data
        standardized_data: list[dict[str, Any]] = []
        for item in data:
            if not isinstance(item, dict):
                continue

            locus_id = item.get("LocusId", "")
            reference_region = item.get("ReferenceRegion", "")
            locus_structure = item.get("LocusStructure", "")

            if not locus_id or not reference_region:
                logger.warning(f"Missing required fields for stripy region: {item}")
                continue

            # Parse reference region (format: "chr:start-end")
            try:
                if ":" not in reference_region or "-" not in reference_region:
                    logger.warning(
                        f"Invalid reference region format: {reference_region}"
                    )
                    continue

                chr_part, coord_part = reference_region.split(":", 1)
                start_str, end_str = coord_part.split("-", 1)

                # Clean chromosome name
                chromosome = chr_part.replace("chr", "")

                # Convert coordinates to integers
                start = int(start_str)
                end = int(end_str)

                # Calculate coverage in base pairs
                coverage = end - start + 1

                # Create comment with locus structure
                comment = f"Repeat structure: {locus_structure}"

                # Default type to "complete" for stripy regions
                region_type = "complete"

                standardized_data.append(
                    {
                        "region_name": locus_id,
                        "chromosome": chromosome,
                        "start": start,
                        "end": end,
                        "coverage": coverage,
                        "source_type": "stripy",
                        "region_type": region_type,
                        "comment": comment,
                    }
                )

            except (ValueError, IndexError) as e:
                logger.warning(f"Error parsing coordinates for {locus_id}: {e}")
                continue

        if not standardized_data:
            logger.warning("No valid stripy regions found")
            return pd.DataFrame()

        result_df = pd.DataFrame(standardized_data)
        logger.info(f"Successfully processed {len(result_df)} stripy repeat regions")
        return result_df

    except Exception as e:
        logger.error(f"Error loading stripy regions: {e}")
        return pd.DataFrame()


def is_available() -> bool:
    """
    Check if regions processing is available.

    Returns:
        True if regions files exist, False otherwise
    """
    manual_file = Path("data/regions/regions_list.xlsx")
    stripy_file = Path("data/regions/stripy_variant_catalog.json")

    return manual_file.exists() or stripy_file.exists()


def get_metadata() -> dict[str, Any]:
    """
    Get metadata about the regions source.

    Returns:
        Dictionary containing source metadata
    """
    return {
        "source_name": "Regions",
        "source_type": "local_files",
        "description": "Manual and stripy repeat genomic regions",
        "data_types": ["manual_regions", "stripy_repeats"],
        "version": "1.0.0",
    }
