"""
PharmGKB Pharmacogenomics SNPs source fetcher.

This module downloads and processes PharmGKB variants data from the PharmGKB API,
filtering for clinically relevant pharmacogenomic variants based on guideline and
clinical annotation counts.

## PharmGKB Data Processing Pipeline

### 1. Data Source
PharmGKB provides comprehensive variant data via their API:
- URL: https://api.pharmgkb.org/v1/download/file/data/variants.zip
- Format: ZIP file containing variants.tsv
- Update frequency: Regularly maintained

### 2. Data Filtering
Variants are filtered based on clinical relevance:
- Guideline Annotation count > 0 (variants with clinical guidelines)
- Level 1/2 Clinical Annotation count > 0 (high-evidence clinical annotations)
- Configurable filter logic (AND/OR)

### 3. Data Structure
Original PharmGKB format:
```
Variant ID    Variant Name    Gene IDs    Gene Symbols    Location
Variant Annotation count    Clinical Annotation count
Level 1/2 Clinical Annotation count    Guideline Annotation count
Label Annotation count    Synonyms
```

### 4. Output Format
Standardized SNP format compatible with pipeline:
```python
{
    'snp': 'rs4244285',           # rsID extracted from Variant Name
    'rsid': 'rs4244285',          # Same as snp for compatibility
    'source': 'PharmGKB',         # Source identifier
    'category': 'pharmacogenomics',# SNP category
    'gene': 'CYP2C19',           # Gene symbol(s)
    'pharmgkb_id': 'PA166156302', # Original PharmGKB Variant ID
    'clinical_annotation_count': 5,
    'guideline_annotation_count': 2,
    'level12_clinical_annotation_count': 3,
    'location': 'NC_000010.11:94842866'
}
```

## Caching Strategy
- Cache directory: `.cache/pharmgkb/`
- Downloaded ZIP: `variants.zip`
- Extracted TSV: `variants.tsv`
- Processed data: `processed_variants.pkl`
- TTL: Configurable (default 7 days)
"""

import logging
import pickle
import re
import zipfile
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import pandas as pd
import requests

logger = logging.getLogger(__name__)


def fetch_pharmacogenomics_snps(config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Fetch pharmacogenomics SNPs from PharmGKB database.

    Downloads the PharmGKB variants.zip file, extracts and filters variants
    based on clinical annotation and guideline annotation counts.

    Returns raw, unharmonized SNP data. Harmonization is handled centrally
    in the Pipeline class for improved efficiency.

    Args:
        config: Configuration dictionary

    Returns:
        DataFrame with raw pharmacogenomics SNPs or None if disabled/failed
    """
    snp_config = config.get("snp_processing", {})

    if not snp_config.get("enabled", False):
        logger.info("SNP processing is disabled")
        return None

    pharmgkb_config = snp_config.get("pharmacogenomics", {})

    if not pharmgkb_config.get("enabled", False):
        logger.info("Pharmacogenomics SNPs are disabled")
        return None

    logger.info("Fetching pharmacogenomics SNPs from PharmGKB")

    try:
        # Get PharmGKB variants data
        variants_df = _download_and_process_pharmgkb_data(pharmgkb_config)

        if variants_df is None or variants_df.empty:
            logger.warning("No PharmGKB variants were successfully processed")
            return None

        # Filter based on clinical relevance
        filtered_df = _filter_clinically_relevant_variants(variants_df, pharmgkb_config)

        if filtered_df.empty:
            logger.warning("No clinically relevant variants found after filtering")
            return None

        # Transform to standard SNP format
        snp_df = _transform_to_snp_format(filtered_df)

        if snp_df.empty:
            logger.warning("No valid SNPs extracted from PharmGKB data")
            return None

        logger.info(
            f"Successfully fetched {len(snp_df)} pharmacogenomics SNPs from PharmGKB"
        )
        return snp_df

    except Exception as e:
        logger.error(f"Failed to fetch PharmGKB pharmacogenomics SNPs: {e}")
        return None


def _download_and_process_pharmgkb_data(config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Download and process PharmGKB variants data with caching.

    Args:
        config: PharmGKB configuration dictionary

    Returns:
        DataFrame with PharmGKB variants or None if failed
    """
    cache_dir = Path(config.get("cache_dir", ".cache/pharmgkb"))
    cache_ttl_days = config.get("cache_ttl_days", 7)
    download_url = config.get(
        "download_url", "https://api.pharmgkb.org/v1/download/file/data/variants.zip"
    )

    # Create cache directory
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Check cache first
    processed_cache = cache_dir / "processed_variants.pkl"
    if _is_cache_valid(processed_cache, cache_ttl_days):
        logger.info("Using cached PharmGKB variants data")
        try:
            with open(processed_cache, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            logger.warning(f"Failed to load cached data: {e}")

    # Download fresh data
    logger.info("Downloading PharmGKB variants data")
    variants_df = _download_pharmgkb_variants(download_url, cache_dir, config)

    if variants_df is not None:
        # Cache the processed data
        try:
            with open(processed_cache, "wb") as f:
                pickle.dump(variants_df, f)
            logger.info("Cached PharmGKB variants data for future use")
        except Exception as e:
            logger.warning(f"Failed to cache processed data: {e}")

    return variants_df


def _download_pharmgkb_variants(
    url: str, cache_dir: Path, config: dict[str, Any]
) -> pd.DataFrame | None:
    """
    Download PharmGKB variants ZIP file and extract TSV data.

    Args:
        url: PharmGKB download URL
        cache_dir: Cache directory path
        config: Configuration dictionary

    Returns:
        DataFrame with variants data or None if failed
    """
    zip_path = cache_dir / "variants.zip"
    cache_dir / "variants.tsv"
    timeout = config.get("timeout", 300)  # 5 minutes default

    try:
        # Download ZIP file
        logger.info(f"Downloading PharmGKB data from {url}")
        response = requests.get(url, timeout=timeout, stream=True)
        response.raise_for_status()

        # Save ZIP file with progress logging
        total_size = int(response.headers.get("content-length", 0))
        downloaded_size = 0

        with open(zip_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded_size += len(chunk)

                    # Log progress every 10MB
                    if (
                        downloaded_size > 0
                        and downloaded_size % (10 * 1024 * 1024) < 8192
                    ):
                        if total_size > 0:
                            progress = (downloaded_size / total_size) * 100
                            logger.info(
                                f"Downloaded {downloaded_size // (1024 * 1024)}MB / "
                                f"{total_size // (1024 * 1024)}MB ({progress:.1f}%)"
                            )
                        else:
                            logger.info(
                                f"Downloaded {downloaded_size // (1024 * 1024)}MB"
                            )

        logger.info(
            f"Successfully downloaded {downloaded_size // (1024 * 1024)}MB ZIP file"
        )

        # Extract TSV file
        logger.info("Extracting variants.tsv from ZIP file")
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            # Look for variants.tsv file
            tsv_files = [
                name for name in zip_ref.namelist() if name.endswith("variants.tsv")
            ]

            if not tsv_files:
                raise ValueError("No variants.tsv file found in ZIP archive")

            tsv_filename = tsv_files[0]
            logger.info(f"Extracting {tsv_filename}")

            with zip_ref.open(tsv_filename) as tsv_file:
                # Read TSV data directly from ZIP
                df = pd.read_csv(tsv_file, sep="\t", low_memory=False)

        logger.info(f"Successfully loaded {len(df)} variants from PharmGKB")
        return df

    except requests.RequestException as e:
        logger.error(f"Failed to download PharmGKB data: {e}")
        return None
    except zipfile.BadZipFile as e:
        logger.error(f"Downloaded file is not a valid ZIP archive: {e}")
        # Clean up corrupted file
        if zip_path.exists():
            zip_path.unlink()
        return None
    except Exception as e:
        logger.error(f"Failed to process PharmGKB data: {e}")
        return None


def _filter_clinically_relevant_variants(
    df: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """
    Filter variants based on clinical relevance criteria.

    Args:
        df: DataFrame with PharmGKB variants
        config: Configuration dictionary with filter settings

    Returns:
        Filtered DataFrame with clinically relevant variants
    """
    filters = config.get("filters", {})
    min_guideline = filters.get("min_guideline_annotations", 1)
    min_level12 = filters.get("min_level12_clinical_annotations", 1)
    filter_logic = filters.get("filter_logic", "OR").upper()

    initial_count = len(df)
    logger.info(f"Filtering {initial_count} variants for clinical relevance")

    # Ensure annotation count columns exist and are numeric
    required_columns = [
        "Guideline Annotation count",
        "Level 1/2 Clinical Annotation count",
    ]

    for col in required_columns:
        if col not in df.columns:
            logger.error(f"Required column '{col}' not found in PharmGKB data")
            return pd.DataFrame()

        # Convert to numeric, replacing non-numeric values with 0
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    # Apply filters
    guideline_mask = df["Guideline Annotation count"] >= min_guideline
    level12_mask = df["Level 1/2 Clinical Annotation count"] >= min_level12

    if filter_logic == "AND":
        final_mask = guideline_mask & level12_mask
        logger.info(
            f"Using AND logic: guideline >= {min_guideline} AND level1/2 >= {min_level12}"
        )
    else:  # OR logic (default)
        final_mask = guideline_mask | level12_mask
        logger.info(
            f"Using OR logic: guideline >= {min_guideline} OR level1/2 >= {min_level12}"
        )

    filtered_df = df[final_mask].copy()
    filtered_count = len(filtered_df)

    logger.info(
        f"Filtered to {filtered_count} clinically relevant variants "
        f"({filtered_count/initial_count*100:.1f}% of total)"
    )

    # Log filter breakdown
    guideline_only = len(df[guideline_mask & ~level12_mask])
    level12_only = len(df[level12_mask & ~guideline_mask])
    both = len(df[guideline_mask & level12_mask])

    logger.info(
        f"Filter breakdown: {guideline_only} guideline-only, "
        f"{level12_only} level1/2-only, {both} both"
    )

    return filtered_df


def _transform_to_snp_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Transform PharmGKB variants to standard SNP format.

    Args:
        df: DataFrame with filtered PharmGKB variants

    Returns:
        DataFrame in standard SNP format
    """
    logger.info("Transforming PharmGKB variants to standard SNP format")

    snp_records = []

    for _, row in df.iterrows():
        try:
            # Extract rsID from Variant Name
            rsid = _extract_rsid(row.get("Variant Name", ""))

            if not rsid:
                continue  # Skip variants without valid rsIDs

            # Create SNP record
            snp_record = {
                "snp": rsid,
                "rsid": rsid,
                "source": "PharmGKB",
                "category": "pharmacogenomics",
                "gene": str(row.get("Gene Symbols", "")).strip(),
                "pharmgkb_id": str(row.get("Variant ID", "")).strip(),
                "clinical_annotation_count": int(
                    row.get("Clinical Annotation count", 0)
                ),
                "guideline_annotation_count": int(
                    row.get("Guideline Annotation count", 0)
                ),
                "level12_clinical_annotation_count": int(
                    row.get("Level 1/2 Clinical Annotation count", 0)
                ),
                "location": str(row.get("Location", "")).strip(),
            }

            snp_records.append(snp_record)

        except Exception as e:
            logger.warning(
                f"Failed to process variant {row.get('Variant ID', 'unknown')}: {e}"
            )
            continue

    if not snp_records:
        logger.warning("No valid SNP records created from PharmGKB data")
        return pd.DataFrame()

    snp_df = pd.DataFrame(snp_records)

    # Remove duplicates based on rsID, keeping the one with highest clinical annotation count
    initial_count = len(snp_df)
    snp_df = snp_df.sort_values("clinical_annotation_count", ascending=False)
    snp_df = snp_df.drop_duplicates(subset=["rsid"], keep="first")
    final_count = len(snp_df)

    if initial_count != final_count:
        logger.info(f"Removed {initial_count - final_count} duplicate rsIDs")

    # Sort by rsID for consistent output
    snp_df = snp_df.sort_values("rsid").reset_index(drop=True)

    logger.info(f"Successfully created {len(snp_df)} SNP records")
    return snp_df


def _extract_rsid(variant_name: str) -> str | None:
    """
    Extract rsID from PharmGKB Variant Name field.

    PharmGKB Variant Name can contain multiple formats:
    - rs12345
    - rs12345; rs67890
    - Other identifiers mixed with rsIDs

    Args:
        variant_name: PharmGKB Variant Name string

    Returns:
        First valid rsID found, or None if no rsID found
    """
    if not variant_name or pd.isna(variant_name):
        return None

    variant_name = str(variant_name).strip()

    # Look for rsIDs using regex
    rsid_pattern = r"\b(rs\d+)\b"
    matches = re.findall(rsid_pattern, variant_name, re.IGNORECASE)

    if matches:
        # Return the first rsID found
        return matches[0].lower()

    return None


def _is_cache_valid(cache_path: Path, ttl_days: int) -> bool:
    """
    Check if cached file is still valid based on TTL.

    Args:
        cache_path: Path to cached file
        ttl_days: Time-to-live in days

    Returns:
        True if cache is valid, False otherwise
    """
    if not cache_path.exists():
        return False

    try:
        # Get file modification time
        mtime = datetime.fromtimestamp(cache_path.stat().st_mtime)

        # Check if file is within TTL
        expiry_time = mtime + timedelta(days=ttl_days)

        return datetime.now() < expiry_time

    except Exception:
        return False


def get_pharmacogenomics_snps_summary(df: pd.DataFrame) -> dict[str, Any]:
    """
    Generate summary statistics for pharmacogenomics SNPs.

    Args:
        df: DataFrame with pharmacogenomics SNPs

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

    # Gene breakdown
    if "gene" in df.columns:
        # Count genes (handle multiple genes per variant)
        all_genes = []
        for genes_str in df["gene"].dropna():
            if genes_str and genes_str != "":
                # Split on common separators
                genes = re.split(r"[;,|]", str(genes_str))
                all_genes.extend([g.strip() for g in genes if g.strip()])

        if all_genes:
            gene_counts = pd.Series(all_genes).value_counts()
            summary["gene_breakdown"] = gene_counts.head(10).to_dict()  # Top 10 genes
            summary["unique_genes"] = len(gene_counts)

    # Clinical annotation statistics
    annotation_cols = [
        "clinical_annotation_count",
        "guideline_annotation_count",
        "level12_clinical_annotation_count",
    ]

    for col in annotation_cols:
        if col in df.columns:
            summary[f"{col}_stats"] = {
                "mean": float(df[col].mean()),
                "max": int(df[col].max()),
                "min": int(df[col].min()),
                "total": int(df[col].sum()),
            }

    # High-evidence variants
    if "guideline_annotation_count" in df.columns:
        summary["with_guidelines"] = int((df["guideline_annotation_count"] > 0).sum())

    if "level12_clinical_annotation_count" in df.columns:
        summary["with_level12_clinical"] = int(
            (df["level12_clinical_annotation_count"] > 0).sum()
        )

    return summary
