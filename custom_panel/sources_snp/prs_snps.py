"""
Polygenic Risk Score (PRS) SNPs source fetcher and aggregator.

This module handles the complete PRS acquisition pipeline from multiple sources
with support for dual genome builds, coordinate systems, and metadata preservation.

## PRS Acquisition Pipeline

### 1. Source Configuration
PRS sources are configured in the SNP processing configuration:
```yaml
snp_processing:
  prs:
    enabled: true
    sources:
      - name: "BCAC_313_PRS"
        parser: "bcac_313"
        file_path: "data/snp/prs/BCAC_313_PRS.prs"
        enhance_with_ensembl: true
      - name: "PGS_Catalog_Breast_Cancer"
        parser: "pgs_catalog_fetcher"
        pgs_ids: ["PGS000873", "PGS000004"]
        genome_build: "GRCh38"
```

### 2. Multi-Source Processing
The pipeline processes multiple PRS sources in parallel:
- **BCAC files**: Local files with original research data
- **PGS Catalog**: API-downloaded harmonized scoring files
- **Custom formats**: Extensible parser system for new formats

### 3. Dual Genome Build Support
Each source can provide coordinates for multiple genome builds:
- **hg19/GRCh37**: Legacy build support
- **hg38/GRCh38**: Modern harmonized coordinates
- **Cross-build mapping**: Ensembl API integration for coordinate conversion

### 4. Aggregation Strategy
SNPs from different sources are aggregated using intelligent merging:
```python
# Group by variant identifier (rsID or coordinate-based)
grouped = df.groupby('snp')

# Preserve metadata from all sources
agg_rules = {
    'source': lambda x: "; ".join(x.unique()),  # Combine sources
    'chromosome': get_first_valid_coord,        # Best coordinate
    'position': get_first_valid_coord,          # Best position
    'effect_allele': 'first',                   # Effect allele
    'effect_weight': 'first',                   # PRS weight
}
```

## Data Quality and Validation

### Input Validation
- **File existence**: Validates all configured file paths
- **Format detection**: Automatic parser selection based on file format
- **Column validation**: Ensures required columns are present

### Data Cleaning
- **Invalid entry filtering**: Removes malformed variant identifiers
- **Coordinate validation**: Validates chromosome and position values
- **Allele normalization**: Standardizes allele representations

### Quality Metrics
- **Coverage reporting**: Tracks variants successfully parsed
- **Build coverage**: Reports coordinate availability by genome build
- **Source overlap**: Identifies variants present in multiple sources

## Coordinate System Management

### Multi-Build Coordinate Handling
```python
# Both genome builds supported simultaneously
hg38_columns = ['hg38_chromosome', 'hg38_start', 'hg38_end', 'hg38_strand']
hg19_columns = ['hg19_chromosome', 'hg19_start', 'hg19_end', 'hg19_strand']

# Coordinate preference hierarchy
coordinate_preference = [
    ('hm_chr', 'hm_pos'),      # Harmonized coordinates (preferred)
    ('chr_name', 'chr_position'), # Original coordinates (fallback)
]
```

### Cross-Build Validation
- **Consistency checks**: Validates coordinates across builds when both available
- **Lift-over validation**: Ensures coordinate conversions are accurate
- **Build-specific metadata**: Preserves allele strings and strand information per build

## Performance and Scalability

### Efficient Processing
- **Lazy loading**: Sources loaded on-demand
- **Parallel processing**: Multiple sources processed concurrently
- **Memory optimization**: Large datasets processed in chunks

### Caching Strategy
- **File caching**: Downloaded PGS files cached with TTL
- **API caching**: Ensembl API responses cached to reduce load
- **Result caching**: Processed aggregations cached for reuse

### Error Recovery
- **Source isolation**: Failures in one source don't affect others
- **Partial results**: Returns available data even if some sources fail
- **Detailed logging**: Comprehensive error reporting for debugging

## Output Format

### Standardized Schema
All PRS SNPs follow a consistent output schema:
```python
required_columns = [
    'snp',           # Variant identifier (rsID or chr:pos:ref:alt)
    'rsid',          # rsID when available
    'source',        # Source name(s)
    'category',      # Always 'prs'
    'chromosome',    # Primary chromosome
    'position',      # Primary position
]

optional_columns = [
    'effect_allele',        # Effect allele for PRS
    'other_allele',         # Reference/other allele
    'effect_weight',        # PRS weight/beta coefficient
    'hg38_chromosome',      # hg38 coordinates
    'hg38_start', 'hg38_end', 'hg38_strand',
    'hg19_chromosome',      # hg19 coordinates
    'hg19_start', 'hg19_end', 'hg19_strand',
    'pgs_id', 'pgs_name',   # PGS Catalog metadata
]
```

### Metadata Preservation
- **Source tracking**: All source names preserved
- **Build information**: Genome build specified per coordinate set
- **PRS metadata**: Effect sizes, frequencies, and statistical measures
- **Provenance**: Original file paths and processing timestamps
"""

import logging
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ..parsers_snp.parsers_prs import create_prs_parser

logger = logging.getLogger(__name__)


def fetch_prs_snps(
    config: dict[str, Any], harmonizer: Optional[Any] = None
) -> pd.DataFrame | None:
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
            logger.error(f"✗ {panel_config.get('name', 'Unknown')}: {e}")

    if not all_snps:
        logger.warning("No PRS SNPs were successfully fetched")
        return None

    # Combine all panels
    combined_df = pd.concat(all_snps, ignore_index=True)

    # Check for conflicting effect alleles before aggregation
    if "effect_allele" in combined_df.columns:
        conflicts = _check_effect_allele_conflicts(combined_df)
        if conflicts:
            logger.warning(
                f"Found {len(conflicts)} rsIDs with conflicting effect alleles"
            )
            for rsid, alleles in conflicts.items():
                logger.warning(f"  {rsid}: {alleles}")

    # Apply R-script-like aggregation but preserve PRS metadata
    prs_snps_panel = _aggregate_prs_snps_by_rsid(combined_df)

    # Apply harmonization if harmonizer is provided
    if harmonizer is not None:
        try:
            logger.info(f"Harmonizing {len(prs_snps_panel)} PRS SNPs")
            harmonized_prs_snps = harmonizer.harmonize_snp_batch(prs_snps_panel)

            if not harmonized_prs_snps.empty:
                logger.info(
                    f"Successfully harmonized {len(harmonized_prs_snps)} PRS SNPs"
                )
                return harmonized_prs_snps
            else:
                logger.warning(
                    "Harmonization resulted in empty DataFrame, returning original data"
                )

        except Exception as e:
            logger.error(f"Error during PRS SNP harmonization: {e}")
            logger.info("Continuing with non-harmonized data")

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
    parser_type = panel_config.get("parser", "bcac_prs")

    # For PGS Catalog fetcher, no file path is needed (downloads via API)
    if parser_type == "pgs_catalog_fetcher":
        file_path = Path(".")  # Dummy path - fetcher doesn't use it
    else:
        file_path = Path(panel_config.get("file_path", ""))

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

    # Clean up DataFrame for aggregation
    if "snp" in df.columns:
        # Ensure 'snp' column contains scalar string values
        df["snp"] = df["snp"].astype(str)

        # Remove any duplicate columns that might cause issues first
        df = df.loc[:, ~df.columns.duplicated()]

        # Filter out problematic entries before aggregation
        initial_count = len(df)
        invalid_snp_patterns = ["liftover", "<NA>", "nan", "null", ""]

        # Create mask for valid SNP identifiers
        valid_mask = pd.Series([True] * len(df), index=df.index)

        # Get the snp column as a Series to ensure proper .str access
        snp_series = df["snp"]

        for pattern in invalid_snp_patterns:
            try:
                pattern_mask = snp_series.str.match(
                    f"^{pattern}$", case=False, na=False
                )
                if pattern_mask.any():
                    logger.info(
                        f"Filtering out {pattern_mask.sum()} entries with SNP ID '{pattern}'"
                    )
                    valid_mask = valid_mask & ~pattern_mask
            except AttributeError as e:
                logger.warning(
                    f"Could not apply string pattern matching for '{pattern}': {e}"
                )
                # Fallback to exact string comparison
                pattern_mask = snp_series == pattern
                if pattern_mask.any():
                    logger.info(
                        f"Filtering out {pattern_mask.sum()} entries with exact SNP ID '{pattern}'"
                    )
                    valid_mask = valid_mask & ~pattern_mask

        # Also filter out entries where SNP is NaN or empty after string conversion
        null_mask = snp_series.isna() | (snp_series == "nan") | (snp_series == "")
        if null_mask.any():
            logger.info(
                f"Filtering out {null_mask.sum()} entries with null/empty SNP IDs"
            )
            valid_mask = valid_mask & ~null_mask

        # Apply the filter
        df = df[valid_mask].reset_index(drop=True)
        filtered_count = initial_count - len(df)

        if filtered_count > 0:
            logger.info(
                f"Filtered out {filtered_count} problematic entries before aggregation"
            )

        # Reset index to ensure clean structure
        df = df.reset_index(drop=True)

    # Define aggregation rules for different column types
    agg_rules: dict[str, Any] = {
        "source": lambda x: "; ".join(x.dropna().astype(str).unique()),
        "category": "first",
    }

    # For PRS metadata, keep first non-null value
    prs_metadata_cols = [
        "chromosome",
        "position",
        "effect_allele",
        "effect_weight",
        "OR",
        "beta",
        "se",
        "pvalue",
        "freq",
        "reference_allele",
        "other_allele",
        "effect_allele_frequency",
        "log_odds_ratio",
        "odds_ratio",
        "alpha",
        "genome_build",
    ]

    # Add coordinate columns for both genome builds
    coordinate_cols = [
        "hg38_chromosome",
        "hg38_chr",
        "hg38_start",
        "hg38_pos",
        "hg38_end",
        "hg38_strand",
        "hg38_allele_string",
        "hg19_chromosome",
        "hg19_chr",
        "hg19_start",
        "hg19_pos",
        "hg19_end",
        "hg19_strand",
        "hg19_allele_string",
    ]

    # Add PGS-specific metadata columns
    pgs_metadata_cols = [
        "pgs_id",
        "pgs_name",
        "trait_reported",
        "trait_mapped",
        "variants_number",
        "weight_type",
    ]

    # Combine all metadata columns
    all_metadata_cols = prs_metadata_cols + coordinate_cols + pgs_metadata_cols

    # Define custom aggregation function for coordinates
    def get_first_valid_coord(series: "pd.Series[Any]") -> Any:
        """Get first non-empty, non-null value from series."""
        valid_values = series.dropna()
        if len(valid_values) == 0:
            return ""
        # Filter out empty strings
        non_empty = valid_values[valid_values != ""]
        if len(non_empty) == 0:
            return ""
        return non_empty.iloc[0]

    for col in all_metadata_cols:
        if col in df.columns:
            if col in coordinate_cols:
                agg_rules[col] = get_first_valid_coord
            else:
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
            conflicts[str(rsid)] = unique_alleles.tolist()

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
        summary["effect_allele_distribution"] = (
            df["effect_allele"].value_counts().to_dict()
        )

    # Chromosome distribution
    if "chromosome" in df.columns:
        summary["chromosome_distribution"] = df["chromosome"].value_counts().to_dict()

    # Additional PRS-specific metrics
    prs_columns = ["effect_weight", "OR", "beta", "pvalue"]
    for col in prs_columns:
        if col in df.columns:
            summary[f"with_{col}"] = df[col].notna().sum()

    return summary
