"""
COSMIC Cancer Gene Census data source extractor.

This module fetches and processes the COSMIC Cancer Gene Census, providing
both germline and somatic evidence scoring based on tier classifications.
"""

import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)


def _download_cosmic_census(url: str, cache_path: Path) -> None:
    """
    Download COSMIC Cancer Gene Census file.

    Args:
        url: URL to download from
        cache_path: Local path to save the file

    Raises:
        requests.RequestException: If download fails
    """
    logger.info(f"Downloading COSMIC Cancer Gene Census from: {url}")

    # Create cache directory
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    # Download with proper headers
    headers = {
        "User-Agent": "custom-panel/1.0 (Python scientific tool for gene panel curation)"
    }

    try:
        response = requests.get(url, headers=headers, timeout=60, stream=True)
        response.raise_for_status()

        # Save file
        with open(cache_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        logger.info(f"Successfully downloaded COSMIC census to: {cache_path}")

    except requests.RequestException as e:
        logger.error(f"Failed to download COSMIC census: {e}")
        raise


def _load_cosmic_census(cache_path: Path) -> pd.DataFrame:
    """
    Load and validate COSMIC Cancer Gene Census file.

    Args:
        cache_path: Path to the cached CSV file

    Returns:
        Loaded DataFrame

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    if not cache_path.exists():
        raise FileNotFoundError(f"COSMIC census file not found: {cache_path}")

    try:
        df = pd.read_csv(cache_path)
        logger.info(f"Loaded COSMIC census with {len(df)} genes")

        # Validate required columns
        required_columns = ["Gene Symbol", "Tier"]
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(
                f"COSMIC census missing required columns: {missing_columns}"
            )

        # Check for germline and somatic columns (may vary by COSMIC version)
        # Try different possible column names
        germline_cols = ["Germline", "Germline Mutation"]
        somatic_cols = ["Somatic", "Somatic Mutation", "Tumour Types(Somatic)"]

        germline_col = None
        somatic_col = None

        for col in germline_cols:
            if col in df.columns:
                germline_col = col
                break

        for col in somatic_cols:
            if col in df.columns:
                somatic_col = col
                break

        if not germline_col and not somatic_col:
            logger.warning("No germline or somatic columns found in COSMIC census")
        else:
            logger.info(
                f"Found germline column: {germline_col}, somatic column: {somatic_col}"
            )

        # Store column names for later use
        df.attrs["germline_col"] = germline_col
        df.attrs["somatic_col"] = somatic_col

        return df

    except Exception as e:
        logger.error(f"Failed to load COSMIC census: {e}")
        raise ValueError(f"Invalid COSMIC census format: {e}") from e


def _is_cache_valid(cache_path: Path, expiry_days: int) -> bool:
    """
    Check if cached file is still valid.

    Args:
        cache_path: Path to cached file
        expiry_days: Number of days before cache expires

    Returns:
        True if cache is valid, False otherwise
    """
    if not cache_path.exists():
        return False

    file_age = datetime.now() - datetime.fromtimestamp(cache_path.stat().st_mtime)
    return file_age < timedelta(days=expiry_days)


def _calculate_cosmic_score(tier: str, tier_weights: dict[str, float]) -> float:
    """
    Calculate evidence score based on COSMIC tier.

    Args:
        tier: COSMIC tier classification
        tier_weights: Mapping of tiers to weights

    Returns:
        Evidence score (0.0-1.0)
    """
    # Clean tier string
    tier = str(tier).strip() if pd.notna(tier) else ""

    # Get weight, defaulting to unknown tier weight
    return tier_weights.get(tier, tier_weights.get("", 0.4))


def _process_cosmic_genes(
    df: pd.DataFrame, category: str, category_config: dict[str, Any]
) -> pd.DataFrame:
    """
    Process COSMIC genes for a specific category (germline or somatic).

    Args:
        df: COSMIC census DataFrame
        category: Either "germline" or "somatic"
        category_config: Configuration for this category

    Returns:
        Standardized DataFrame for this category
    """
    if not category_config.get("enabled", False):
        logger.info(f"COSMIC {category} scoring is disabled")
        return pd.DataFrame()

    # Get the appropriate column name
    col_name = df.attrs.get(f"{category}_col")
    if not col_name:
        logger.warning(f"No {category} column found in COSMIC census")
        return pd.DataFrame()

    # Filter for genes relevant to this category
    # Look for "yes", "y", or non-empty values indicating presence
    category_df = df[df[col_name].notna()].copy()
    category_df = category_df[
        category_df[col_name].astype(str).str.lower().isin(["yes", "y"])
        | (category_df[col_name].astype(str).str.len() > 0)
    ]

    if category_df.empty:
        logger.warning(f"No {category} genes found in COSMIC census")
        return pd.DataFrame()

    logger.info(f"Found {len(category_df)} COSMIC {category} genes")

    # Calculate evidence scores based on tier weights
    tier_weights = category_config.get("tier_weights", {"": 0.4})
    evidence_scores = [
        _calculate_cosmic_score(tier, tier_weights) for tier in category_df["Tier"]
    ]

    # Create source details
    source_details = [
        f"Tier:{tier}|Category:{category}|Date:{datetime.now().strftime('%Y-%m-%d')}"
        for tier in category_df["Tier"]
    ]

    # Create standardized DataFrame
    genes = category_df["Gene Symbol"].tolist()
    source_name = f"COSMIC_{category.title()}"

    result_df = create_standard_dataframe(
        genes=genes,
        source_name=source_name,
        evidence_scores=evidence_scores,
        source_details=source_details,
        gene_names_reported=genes,
    )

    return result_df


def fetch_cosmic_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch COSMIC Cancer Gene Census data with caching and dual-category processing.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with COSMIC data for both germline and somatic categories
    """
    cosmic_config = config.get("data_sources", {}).get("cosmic", {})

    if not cosmic_config.get("enabled", False):
        logger.info("COSMIC data source is disabled")
        return pd.DataFrame()

    # Get configuration parameters
    census_url = cosmic_config.get("census_url")
    if not census_url:
        logger.error("COSMIC census_url not configured")
        return pd.DataFrame()

    cache_dir = Path(cosmic_config.get("cache_dir", ".cache/cosmic"))
    cache_expiry_days = cosmic_config.get("cache_expiry_days", 30)
    cache_path = cache_dir / "cosmic_gene_census.csv"

    # Check cache and download if needed
    if not _is_cache_valid(cache_path, cache_expiry_days):
        try:
            _download_cosmic_census(census_url, cache_path)
        except requests.RequestException as e:
            logger.error(f"Failed to download COSMIC census: {e}")
            if cache_path.exists():
                logger.warning("Using expired cache file")
            else:
                logger.error("No cache file available, cannot proceed")
                return pd.DataFrame()
    else:
        logger.info(f"Using cached COSMIC census: {cache_path}")

    # Load the census data
    try:
        df = _load_cosmic_census(cache_path)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Failed to load COSMIC census: {e}")
        return pd.DataFrame()

    # Process both germline and somatic categories
    result_dfs = []

    # Process germline
    germline_config = cosmic_config.get("germline_scoring", {})
    germline_df = _process_cosmic_genes(df, "germline", germline_config)
    if not germline_df.empty:
        result_dfs.append(germline_df)

    # Process somatic
    somatic_config = cosmic_config.get("somatic_scoring", {})
    somatic_df = _process_cosmic_genes(df, "somatic", somatic_config)
    if not somatic_df.empty:
        result_dfs.append(somatic_df)

    # Combine results
    if not result_dfs:
        logger.warning("No COSMIC data processed for any category")
        return pd.DataFrame()

    final_df = pd.concat(result_dfs, ignore_index=True)
    logger.info(f"Created COSMIC dataset with {len(final_df)} total gene records")

    return final_df


def validate_cosmic_config(config: dict[str, Any]) -> list[str]:
    """
    Validate COSMIC configuration.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation errors
    """
    errors: list[str] = []
    cosmic_config = config.get("data_sources", {}).get("cosmic", {})

    if not cosmic_config.get("enabled", False):
        return errors  # Skip validation if disabled

    # Check required parameters
    if not cosmic_config.get("census_url"):
        errors.append("COSMIC census_url is required")

    # Validate cache settings
    cache_expiry = cosmic_config.get("cache_expiry_days")
    if cache_expiry is not None and (
        not isinstance(cache_expiry, int) or cache_expiry < 1
    ):
        errors.append("COSMIC cache_expiry_days must be a positive integer")

    # Validate scoring configurations
    for category in ["germline_scoring", "somatic_scoring"]:
        category_config = cosmic_config.get(category, {})
        if category_config.get("enabled", False):
            tier_weights = category_config.get("tier_weights", {})
            if not isinstance(tier_weights, dict):
                errors.append(f"COSMIC {category} tier_weights must be a dictionary")
            else:
                for tier, weight in tier_weights.items():
                    if not isinstance(weight, int | float) or weight < 0 or weight > 1:
                        errors.append(
                            f"COSMIC {category} tier weight for '{tier}' must be between 0 and 1"
                        )

    return errors


def get_cosmic_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary of COSMIC configuration and cached data.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    cosmic_config = config.get("data_sources", {}).get("cosmic", {})

    summary = {
        "enabled": cosmic_config.get("enabled", False),
        "census_url": cosmic_config.get("census_url"),
        "cache_dir": cosmic_config.get("cache_dir", ".cache/cosmic"),
        "cache_expiry_days": cosmic_config.get("cache_expiry_days", 30),
        "germline_enabled": cosmic_config.get("germline_scoring", {}).get(
            "enabled", False
        ),
        "somatic_enabled": cosmic_config.get("somatic_scoring", {}).get(
            "enabled", False
        ),
        "validation_errors": validate_cosmic_config(config),
    }

    # Check cache status
    cache_path = Path(summary["cache_dir"]) / "cosmic_gene_census.csv"
    summary["cache_exists"] = cache_path.exists()
    if cache_path.exists():
        cache_age = datetime.now() - datetime.fromtimestamp(cache_path.stat().st_mtime)
        summary["cache_age_days"] = cache_age.days
        summary["cache_valid"] = cache_age < timedelta(
            days=summary["cache_expiry_days"]
        )
    else:
        summary["cache_age_days"] = None
        summary["cache_valid"] = False

    return summary
