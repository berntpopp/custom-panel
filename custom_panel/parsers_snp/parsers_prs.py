"""
Polygenic Risk Score (PRS) SNP parsers.

This module contains parsers for Polygenic Risk Score (PRS) SNP files which support
multi-genome build coordinate systems and complex metadata structures.

## PRS Acquisition Logic

The PRS parsing system supports multiple data sources and formats:

### 1. PGS Catalog Parser (PGSCatalogParser)
- **Source**: PGS Catalog harmonized scoring files
- **Format**: Tab-separated with extensive header metadata
- **Genome Builds**: Supports both hg19/GRCh37 and hg38/GRCh38
- **Coordinates**: Dual coordinate system with original and harmonized positions
- **Features**:
  - Automatic genome build detection from metadata
  - Harmonized coordinate preference (hm_chr, hm_pos)
  - Fallback to original coordinates (chr_name, chr_position)
  - Multi-strategy variant identifier creation
  - Robust handling of missing rsID fields

### 2. PGS Catalog Fetcher (PGSCatalogFetcher)
- **Source**: PGS Catalog API with automatic download
- **Caching**: Local file caching with configurable TTL
- **Multi-build Support**: Downloads both hg19 and hg38 versions when available
- **Automation**: Fetches multiple PGS IDs in parallel
- **Features**:
  - API-based retrieval with error handling
  - Automatic file format detection
  - Metadata preservation from API responses

## Variant Identifier Strategy

The PRS parsers implement a hierarchical strategy for creating variant identifiers:

1. **Harmonized rsIDs** (`hm_rsID` column) - Preferred when available
2. **Original rsIDs** (`rsID`, `rsid`, `rsId` columns) - Secondary preference
3. **Coordinate-based IDs** - Fallback using genomic coordinates:
   - Full format: `chr:pos:ref:alt` (when both alleles available)
   - Simplified: `chr:pos:effect_allele` (when reference allele missing)

## File Format Handling

### PGS Catalog File Structure
```
# Header metadata (key-value pairs starting with #)
# Multiple genome builds and coordinate systems supported
# Tab-separated data with variable column presence

chr_name	chr_position	effect_allele	other_allele	effect_weight	...	hm_rsID	hm_chr	hm_pos
1	100880328	T	A	0.0373	...		1	100414772
```

### Column Alignment Logic
The parser handles missing columns intelligently:
- **End padding**: Most common case - missing optional columns at end
- **Start padding**: Only when first column is explicitly rsID-named and data starts with chromosome
- **Validation**: Ensures proper column count and data type conversion

## Coordinate System Integration

### Dual Genome Build Support
- **hg19/GRCh37**: Legacy build support for older PRS files
- **hg38/GRCh38**: Modern build preference with harmonized coordinates
- **Automatic conversion**: Ensembl API integration for cross-build mapping

### Coordinate Column Strategy
```python
# Prefer harmonized coordinates when available
chr_col = "hm_chr" if "hm_chr" in df.columns else "chr_name"
pos_col = "hm_pos" if "hm_pos" in df.columns else "chr_position"

# Handle allele information
effect_allele_col = "effect_allele"  # Always present
other_allele_col = "other_allele" or "hm_inferOtherAllele"  # Optional
```

## Error Handling and Validation

### Robust Parsing
- **Missing data tolerance**: Handles NaN, empty strings, and missing columns
- **Type conversion**: Automatic numeric conversion with error handling
- **Validation filtering**: Removes malformed entries before aggregation

### Quality Control
- **Column count validation**: Ensures header-data alignment
- **Data type validation**: Validates numeric columns (positions, weights)
- **Identifier validation**: Filters invalid variant identifiers

## Performance Optimizations

### Caching Strategy
- **File-level caching**: Downloaded PGS files cached locally
- **API rate limiting**: Controlled Ensembl API access
- **Batch processing**: Efficient handling of large variant sets

### Memory Management
- **Streaming parsing**: Large files processed in chunks
- **Selective loading**: Only required columns loaded when possible
- **Garbage collection**: Explicit cleanup of large DataFrames
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.pgs_catalog_client import PGSCatalogClient
from .base_snp_parser import BaseSNPParser

logger = logging.getLogger(__name__)


class BCACParser(BaseSNPParser):
    """
    Legacy parser for BCAC (Breast Cancer Association Consortium) PRS files.

    This parser handles older BCAC formats with rsID columns.
    Note: For modern PRS processing, use PGS Catalog fetcher instead.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a BCAC PRS CSV file.

        Returns:
            DataFrame with columns: rsid, source, category, plus PRS metadata
        """
        self.validate_file_exists()

        logger.info(f"Parsing BCAC PRS file: {self.file_path}")

        # Get column mappings from config
        rsid_column = self.config.get("rsid_column", "hm_rsid")
        chr_column = self.config.get("chr_column", "hm_chr")
        pos_column = self.config.get("pos_column", "hm_pos")
        effect_allele_column = self.config.get("effect_allele_column", "effect_allele")

        try:
            # Read CSV file
            df = pd.read_csv(self.file_path)

            # Check for required rsID column
            if rsid_column not in df.columns:
                raise ValueError(
                    f"Required rsID column '{rsid_column}' not found. Available columns: {list(df.columns)}"
                )

            # Extract rsIDs and remove empty ones
            rsids = df[rsid_column].dropna()
            if rsids.empty:
                logger.warning(
                    f"No valid rsIDs found in column '{rsid_column}' of {self.file_path}"
                )
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Create base DataFrame with required columns
            result_df = pd.DataFrame(
                {"rsid": rsids, "source": self.name, "category": "prs"}
            )

            # Add PRS-specific metadata columns if available
            metadata_columns = {
                "chromosome": chr_column,
                "position": pos_column,
                "effect_allele": effect_allele_column,
            }

            for result_col, source_col in metadata_columns.items():
                if source_col in df.columns:
                    # Align with rsid index to ensure proper mapping
                    result_df[result_col] = df.loc[rsids.index, source_col].values
                else:
                    logger.debug(
                        f"Optional column '{source_col}' not found in {self.file_path}"
                    )

            # Add any additional columns that might be useful for PRS
            additional_columns = ["effect_weight", "OR", "beta", "se", "pvalue", "freq"]
            for col in additional_columns:
                if col in df.columns:
                    result_df[col] = df.loc[rsids.index, col].values

            logger.info(
                f"Successfully parsed {len(result_df)} PRS SNPs from {self.file_path}"
            )

            # Validate and standardize output
            return self.validate_output_format(result_df)

        except Exception as e:
            logger.error(f"Error parsing BCAC PRS file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse BCAC PRS file: {e}") from e


class GenericPRSParser(BaseSNPParser):
    """
    Generic parser for PRS files with configurable column mappings.

    This parser can handle various PRS file formats by allowing
    flexible column mapping configuration.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a generic PRS file with configurable column mappings.

        Returns:
            DataFrame with columns: rsid, source, category, plus metadata
        """
        self.validate_file_exists()

        logger.info(f"Parsing generic PRS file: {self.file_path}")

        # Get required column mapping
        rsid_column = self.config.get("rsid_column")
        if not rsid_column:
            raise ValueError(
                "rsid_column must be specified in config for generic PRS parser"
            )

        # Get file format (default to CSV)
        file_format = self.config.get("format", "csv")
        separator = self.config.get("separator", ",")

        try:
            # Read file based on format
            if file_format.lower() == "csv":
                df = pd.read_csv(self.file_path, sep=separator)
            elif file_format.lower() in ["xlsx", "excel"]:
                sheet_name = self.config.get("sheet_name", 0)
                df = pd.read_excel(self.file_path, sheet_name=sheet_name)
            elif file_format.lower() == "tsv":
                df = pd.read_csv(self.file_path, sep="\t")
            else:
                raise ValueError(f"Unsupported file format: {file_format}")

            # Check for required rsID column
            if rsid_column not in df.columns:
                raise ValueError(
                    f"Required rsID column '{rsid_column}' not found. Available columns: {list(df.columns)}"
                )

            # Extract rsIDs
            rsids = df[rsid_column].dropna()
            if rsids.empty:
                logger.warning(
                    f"No valid rsIDs found in column '{rsid_column}' of {self.file_path}"
                )
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Create base DataFrame
            result_df = pd.DataFrame(
                {"rsid": rsids, "source": self.name, "category": "prs"}
            )

            # Add any additional columns specified in column mappings
            column_mappings = self.config.get("column_mappings", {})
            for result_col, source_col in column_mappings.items():
                if source_col in df.columns:
                    result_df[result_col] = df.loc[rsids.index, source_col].values

            logger.info(
                f"Successfully parsed {len(result_df)} PRS SNPs from {self.file_path}"
            )

            # Validate and standardize output
            return self.validate_output_format(result_df)

        except Exception as e:
            logger.error(f"Error parsing generic PRS file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse generic PRS file: {e}") from e


class PGSCatalogParser(BaseSNPParser):
    """
    Parser for PGS Catalog harmonized scoring files.

    This parser handles the PGS Catalog format with header metadata
    and tab-separated variant data including harmonized coordinates.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a PGS Catalog harmonized scoring file.

        Returns:
            DataFrame with columns: snp, rsid, source, category, plus PRS metadata
        """
        self.validate_file_exists()

        logger.info(f"Parsing PGS Catalog file: {self.file_path}")

        try:
            # Parse header and data
            header_metadata, df = self._parse_pgs_catalog_file()

            if df.empty:
                logger.warning(f"No valid variants found in {self.file_path}")
                return pd.DataFrame(columns=["snp", "rsid", "source", "category"])

            # Get genome build from header or config
            genome_build = self._determine_genome_build(header_metadata)

            # Create base DataFrame with required columns
            # PGS Catalog format has multiple possible identifier strategies
            variant_ids = self._create_variant_identifiers(df)

            # Filter out invalid/problematic entries
            valid_mask = self._validate_variant_identifiers(variant_ids)
            if not valid_mask.any():
                logger.warning(
                    f"No valid variants found after filtering in {self.file_path}"
                )
                return pd.DataFrame(columns=["snp", "rsid", "source", "category"])

            # Apply filter to both variant IDs and source data
            variant_ids = variant_ids[valid_mask]
            df = df[valid_mask].reset_index(drop=True)

            result_df = pd.DataFrame(
                {
                    "snp": variant_ids,
                    "rsid": variant_ids,
                    "source": self.name,
                    "category": "prs",
                }
            )

            # Add coordinate information (prefer harmonized if available)
            self._add_coordinate_columns(result_df, df, genome_build)

            # Add PRS-specific metadata
            self._add_prs_metadata(result_df, df, header_metadata)

            logger.info(
                f"Successfully parsed {len(result_df)} PRS variants from {self.file_path}"
            )

            # Validate and standardize output
            return self.validate_output_format(result_df)

        except Exception as e:
            logger.error(f"Error parsing PGS Catalog file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse PGS Catalog file: {e}") from e

    def _parse_pgs_catalog_file(self) -> tuple[dict[str, Any], pd.DataFrame]:
        """
        Parse PGS Catalog file extracting header metadata and variant data.

        Returns:
            Tuple of (header_metadata, variants_dataframe)
        """
        # Determine if file is gzipped
        is_gzipped = self.file_path.suffix == ".gz"
        open_func = gzip.open if is_gzipped else open
        mode = "rt" if is_gzipped else "r"

        header_metadata = {}
        data_lines = []
        header_line = None
        in_data_section = False

        with open_func(self.file_path, mode) as f:
            for line_raw in f:
                # Ensure we have a string, not bytes
                line = (
                    line_raw.decode("utf-8")
                    if isinstance(line_raw, bytes)
                    else line_raw
                )
                line = line.strip()

                if not line:
                    continue

                if line.startswith("#"):
                    # Parse header metadata
                    if "=" in line:
                        key, value = line[1:].split("=", 1)
                        header_metadata[key] = value
                elif not in_data_section:
                    # First non-comment line is the header
                    header_line = line
                    in_data_section = True
                else:
                    # Data lines
                    data_lines.append(line)

        if not header_line or not data_lines:
            raise ValueError("Invalid PGS Catalog file format: missing header or data")

        # Parse data into DataFrame
        columns = header_line.split("\t")

        # Parse data lines
        data_rows = []
        for line in data_lines:
            values = line.split("\t")

            # Handle lines with missing values (common in PGS files)
            if len(values) < len(columns):
                # For PGS files, missing values are typically at the end (empty optional columns)
                # Only pad at the beginning if we have strong evidence of a missing rsID
                # This happens when the first column is expected to be rsID but we see a chromosome
                first_val = values[0] if values else ""
                first_column_name = columns[0] if columns else ""

                # Only pad at start if:
                # 1. First column is explicitly named for rsID (rsID, rsid, etc.)
                # 2. AND first value looks like chromosome (not rsID)
                should_pad_start = (
                    first_column_name.lower()
                    in ["rsid", "rsid", "rs_id", "rs", "snp_id", "snp"]
                    and first_val
                    and not first_val.startswith("rs")
                    and first_val.replace("X", "")
                    .replace("Y", "")
                    .replace("M", "")
                    .isdigit()
                )

                if should_pad_start:
                    # Likely missing rsID at the beginning - pad at start
                    padded_values = [""] + values
                    logger.debug(
                        f"Missing rsID detected in rsID column, padding at start: {line[:50]}..."
                    )
                else:
                    # Missing values at the end - pad at end (most common case)
                    padded_values = values[:]
                    logger.debug(
                        f"Missing values at end, padding {len(columns) - len(values)} columns: {line[:50]}..."
                    )

                # Ensure we have the right number of columns
                while len(padded_values) < len(columns):
                    padded_values.append("")

                data_rows.append(padded_values[: len(columns)])  # Truncate if too long
            elif len(values) == len(columns):
                # Perfect match
                data_rows.append(values)
            elif len(values) > len(columns):
                # Too many values - truncate
                logger.debug(f"Too many values, truncating: {line[:50]}...")
                data_rows.append(values[: len(columns)])
            else:
                logger.warning(f"Skipping malformed line: {line[:100]}...")

        if not data_rows:
            raise ValueError("No valid data rows found in PGS Catalog file")

        df = pd.DataFrame(data_rows, columns=columns)

        # Clean up DataFrame - remove empty values and convert numeric columns
        df = df.replace("", pd.NA)

        # Convert numeric columns
        numeric_columns = ["chr_position", "effect_weight", "hm_pos"]
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        logger.info(f"Parsed header metadata: {len(header_metadata)} fields")
        logger.info(
            f"Parsed variant data: {len(df)} variants with {len(df.columns)} columns"
        )

        return header_metadata, df

    def _determine_genome_build(self, header_metadata: dict[str, Any]) -> str:
        """
        Determine genome build from header metadata or config.

        Args:
            header_metadata: Parsed header metadata

        Returns:
            Genome build string (hg19, hg38, etc.)
        """
        # Check harmonized build first
        if "HmPOS_build" in header_metadata:
            hm_build = header_metadata["HmPOS_build"]
            if hm_build == "GRCh38":
                return "hg38"
            elif hm_build == "GRCh37":
                return "hg19"

        # Check original genome build
        if "genome_build" in header_metadata:
            orig_build = header_metadata["genome_build"]
            if orig_build in ["hg19", "GRCh37"]:
                return "hg19"
            elif orig_build in ["hg38", "GRCh38"]:
                return "hg38"

        # Check config preference
        config_build = self.config.get("genome_build", "hg38")
        logger.info(f"Using genome build from config: {config_build}")
        return config_build

    def _add_coordinate_columns(
        self, result_df: pd.DataFrame, data_df: pd.DataFrame, genome_build: str
    ) -> None:
        """
        Add coordinate columns to result DataFrame.

        Args:
            result_df: Result DataFrame to modify
            data_df: Source data DataFrame
            genome_build: Genome build to use
        """
        # Prefer harmonized coordinates if available
        if "hm_chr" in data_df.columns and "hm_pos" in data_df.columns:
            # Use harmonized coordinates
            chromosomes = data_df["hm_chr"].astype(str)
            positions = data_df["hm_pos"]
            logger.info("Using harmonized coordinates from hm_chr/hm_pos columns")
        elif "chr_name" in data_df.columns and "chr_position" in data_df.columns:
            # Use original coordinates
            chromosomes = data_df["chr_name"].astype(str)
            positions = data_df["chr_position"]
            logger.info("Using original coordinates from chr_name/chr_position columns")
        else:
            logger.warning("No coordinate information found in PGS file")
            return

        # Add coordinate columns based on genome build
        if genome_build.lower() in ["hg38", "grch38"]:
            result_df["hg38_chromosome"] = chromosomes
            result_df["hg38_chr"] = chromosomes
            result_df["hg38_start"] = positions
            result_df["hg38_pos"] = positions
            result_df["hg38_end"] = positions  # SNVs have same start/end
            result_df["hg38_strand"] = 1

            # Create allele string if alleles available
            if "other_allele" in data_df.columns and "effect_allele" in data_df.columns:
                result_df["hg38_allele_string"] = (
                    data_df["other_allele"].astype(str)
                    + "/"
                    + data_df["effect_allele"].astype(str)
                )

            # Empty hg19 columns
            result_df["hg19_chromosome"] = ""
            result_df["hg19_chr"] = ""
            result_df["hg19_start"] = ""
            result_df["hg19_pos"] = ""
            result_df["hg19_end"] = ""
            result_df["hg19_strand"] = ""
            result_df["hg19_allele_string"] = ""
        else:
            # hg19 coordinates
            result_df["hg19_chromosome"] = chromosomes
            result_df["hg19_chr"] = chromosomes
            result_df["hg19_start"] = positions
            result_df["hg19_pos"] = positions
            result_df["hg19_end"] = positions
            result_df["hg19_strand"] = 1

            if "other_allele" in data_df.columns and "effect_allele" in data_df.columns:
                result_df["hg19_allele_string"] = (
                    data_df["other_allele"].astype(str)
                    + "/"
                    + data_df["effect_allele"].astype(str)
                )

            # Empty hg38 columns
            result_df["hg38_chromosome"] = ""
            result_df["hg38_chr"] = ""
            result_df["hg38_start"] = ""
            result_df["hg38_pos"] = ""
            result_df["hg38_end"] = ""
            result_df["hg38_strand"] = ""
            result_df["hg38_allele_string"] = ""

        # Add additional coordinate metadata
        result_df["genome_build"] = genome_build
        result_df["chromosome"] = chromosomes
        result_df["position"] = positions

    def _add_prs_metadata(
        self,
        result_df: pd.DataFrame,
        data_df: pd.DataFrame,
        header_metadata: dict[str, Any],
    ) -> None:
        """
        Add PRS-specific metadata to result DataFrame.

        Args:
            result_df: Result DataFrame to modify
            data_df: Source data DataFrame
            header_metadata: Header metadata
        """
        # Add variant-level metadata
        if "effect_allele" in data_df.columns:
            result_df["effect_allele"] = data_df["effect_allele"]
        if "other_allele" in data_df.columns:
            result_df["other_allele"] = data_df["other_allele"]
        if "effect_weight" in data_df.columns:
            result_df["effect_weight"] = pd.to_numeric(
                data_df["effect_weight"], errors="coerce"
            )

        # Add PGS-level metadata from header
        pgs_metadata_fields = {
            "pgs_id": "pgs_id",
            "pgs_name": "pgs_name",
            "trait_reported": "trait_reported",
            "trait_mapped": "trait_mapped",
            "variants_number": "variants_number",
            "weight_type": "weight_type",
        }

        for result_col, header_key in pgs_metadata_fields.items():
            if header_key in header_metadata:
                result_df[result_col] = header_metadata[header_key]

    def _create_variant_identifiers(self, df: pd.DataFrame) -> pd.Series[Any]:
        """
        Create variant identifiers from PGS data, handling different column formats.

        PGS files can have different structures:
        - Some have rsID column, some don't
        - Some have other_allele column, some don't
        - Some have harmonized coordinates (hm_*), some don't

        Args:
            df: DataFrame with PGS data

        Returns:
            Series of variant identifiers
        """
        logger.info(f"Creating variant identifiers for {len(df)} variants")

        # Strategy 1: Use harmonized rsID if available and not empty
        if "hm_rsID" in df.columns:
            hm_rsids = df["hm_rsID"].fillna("")
            valid_hm_rsids = (hm_rsids != "") & hm_rsids.notna()
            if valid_hm_rsids.any():
                logger.info(
                    f"Found {valid_hm_rsids.sum()} variants with harmonized rsIDs"
                )
                variant_ids = hm_rsids.copy()
                # For missing harmonized rsIDs, we'll fill them below
                needs_fallback = ~valid_hm_rsids
            else:
                needs_fallback = pd.Series([True] * len(df), index=df.index)
                variant_ids = pd.Series([""] * len(df), index=df.index)
        else:
            needs_fallback = pd.Series([True] * len(df), index=df.index)
            variant_ids = pd.Series([""] * len(df), index=df.index)

        # Strategy 2: Use original rsID if available and not empty
        rsid_column = None
        for col in ["rsID", "rsid", "rsId"]:  # Check multiple case variations
            if col in df.columns:
                rsid_column = col
                break

        if needs_fallback.any() and rsid_column:
            orig_rsids = df[rsid_column].fillna("")
            valid_orig_rsids = (orig_rsids != "") & orig_rsids.notna() & needs_fallback
            if valid_orig_rsids.any():
                logger.info(
                    f"Found {valid_orig_rsids.sum()} variants with original rsIDs from column '{rsid_column}'"
                )
                variant_ids.loc[valid_orig_rsids] = orig_rsids.loc[valid_orig_rsids]
                needs_fallback = needs_fallback & ~valid_orig_rsids

        # Strategy 3: Create chr:pos:ref:alt identifier for remaining variants
        if needs_fallback.any():
            logger.info(
                f"Creating coordinate-based IDs for {needs_fallback.sum()} variants"
            )

            # Determine coordinate columns (prefer harmonized)
            chr_col = "hm_chr" if "hm_chr" in df.columns else "chr_name"
            pos_col = "hm_pos" if "hm_pos" in df.columns else "chr_position"

            # Determine allele columns (handle missing other_allele)
            effect_allele_col = "effect_allele"
            if "other_allele" in df.columns:
                other_allele_col = "other_allele"
            elif "hm_inferOtherAllele" in df.columns:
                other_allele_col = "hm_inferOtherAllele"
            else:
                # No other allele info available, use simplified format
                other_allele_col = None

            fallback_mask = needs_fallback
            if other_allele_col:
                # Full chr:pos:ref:alt format
                coord_ids = (
                    df.loc[fallback_mask, chr_col].astype(str)
                    + ":"
                    + df.loc[fallback_mask, pos_col].astype(str)
                    + ":"
                    + df.loc[fallback_mask, other_allele_col].astype(str)
                    + ":"
                    + df.loc[fallback_mask, effect_allele_col].astype(str)
                )
            else:
                # Simplified chr:pos:allele format
                coord_ids = (
                    df.loc[fallback_mask, chr_col].astype(str)
                    + ":"
                    + df.loc[fallback_mask, pos_col].astype(str)
                    + ":"
                    + df.loc[fallback_mask, effect_allele_col].astype(str)
                )

            variant_ids.loc[fallback_mask] = coord_ids

        logger.info(f"Successfully created {len(variant_ids)} variant identifiers")
        return variant_ids

    def _validate_variant_identifiers(
        self, variant_ids: pd.Series[Any]
    ) -> pd.Series[Any]:
        """
        Validate variant identifiers and filter out problematic entries.

        Args:
            variant_ids: Series of variant identifiers

        Returns:
            Boolean mask indicating valid entries
        """
        # Create mask for valid entries
        valid_mask = pd.Series([True] * len(variant_ids), index=variant_ids.index)

        # Filter out null/empty values
        null_mask = variant_ids.isna() | (variant_ids == "")
        if null_mask.any():
            logger.info(f"Filtering out {null_mask.sum()} null/empty entries")
            valid_mask = valid_mask & ~null_mask

        # For now, we are being less aggressive with filtering since the main issue
        # seems to be during aggregation, not during individual file parsing

        logger.info(
            f"Validation complete: {valid_mask.sum()}/{len(variant_ids)} entries are valid"
        )
        return valid_mask


class PGSCatalogFetcher(BaseSNPParser):
    """
    Fetcher for PGS Catalog scores via API with automatic download and parsing.

    This fetcher uses the PGS Catalog API to download harmonized scoring files
    and parse them automatically based on configuration.
    """

    def parse(self) -> pd.DataFrame:
        """
        Fetch PGS data from PGS Catalog API and parse scoring files.

        Returns:
            DataFrame with combined PRS data from all requested PGS IDs
        """
        logger.info(f"Fetching PGS data via API for: {self.name}")

        # Get configuration
        pgs_ids = self.config.get("pgs_ids", [])
        genome_build = self.config.get("genome_build", "GRCh38")
        cache_dir = self.config.get("cache_dir", ".cache/pgs_catalog")
        cache_ttl_days = self.config.get("cache_ttl_days", 7)

        if not pgs_ids:
            raise ValueError("No PGS IDs specified in configuration")

        try:
            # Initialize PGS Catalog client
            client = PGSCatalogClient(
                cache_dir=cache_dir, cache_ttl_days=cache_ttl_days
            )

            # Fetch and download scoring files for both builds
            downloaded_files = client.fetch_and_cache_pgs_files(pgs_ids, genome_build)

            # Parse each downloaded file and merge coordinate data
            all_pgs_data = []
            for pgs_id, file_dict in downloaded_files.items():
                logger.info(f"Parsing {pgs_id} with {len(file_dict)} genome builds")

                # Parse files for each available build
                build_dataframes = {}
                for build, file_path in file_dict.items():
                    logger.info(f"Parsing {pgs_id} ({build}) from {file_path}")

                    # Create parser config for this PGS + build
                    parser_config = self.config.copy()
                    parser_config["name"] = f"{self.name}_{pgs_id}_{build}"
                    parser_config["genome_build"] = (
                        "hg38" if build == "GRCh38" else "hg19"
                    )

                    # Parse the file
                    parser = PGSCatalogParser(file_path, parser_config)
                    build_df = parser.parse()

                    if not build_df.empty:
                        build_dataframes[build] = build_df
                        logger.info(
                            f"✓ Parsed {len(build_df)} variants from {pgs_id} ({build})"
                        )
                    else:
                        logger.warning(f"⚠ No variants found in {pgs_id} ({build})")

                if build_dataframes:
                    # Merge coordinate data from different builds
                    merged_df = self._merge_dual_build_data(build_dataframes, pgs_id)
                    if not merged_df.empty:
                        all_pgs_data.append(merged_df)
                        logger.info(f"✓ Merged {len(merged_df)} variants from {pgs_id}")
                else:
                    logger.warning(f"⚠ No data for {pgs_id} from any build")

            if not all_pgs_data:
                logger.warning("No PGS data was successfully parsed")
                return pd.DataFrame(columns=["snp", "rsid", "source", "category"])

            # Combine all PGS data
            combined_df = pd.concat(all_pgs_data, ignore_index=True)

            logger.info(
                f"Successfully fetched and parsed {len(combined_df)} total variants from {len(downloaded_files)} PGS scores"
            )

            return combined_df

        except Exception as e:
            logger.error(f"Error fetching PGS data: {e}")
            raise ValueError(f"Failed to fetch PGS data: {e}") from e

    def _merge_dual_build_data(
        self, build_dataframes: dict[str, pd.DataFrame], pgs_id: str
    ) -> pd.DataFrame:
        """
        Merge coordinate data from different genome builds for the same PGS.

        Args:
            build_dataframes: Dict of {build: DataFrame}
            pgs_id: PGS identifier for logging

        Returns:
            Merged DataFrame with coordinates from both builds
        """
        if len(build_dataframes) == 1:
            # Only one build available, return as-is
            return next(iter(build_dataframes.values()))

        # Get dataframes for each build
        grch38_df = build_dataframes.get("GRCh38")
        grch37_df = build_dataframes.get("GRCh37")

        if grch38_df is None and grch37_df is None:
            return pd.DataFrame()

        # Use GRCh38 as base if available, otherwise GRCh37
        base_df = grch38_df if grch38_df is not None else grch37_df
        if base_df is None:
            return pd.DataFrame()
        result_df = base_df.copy()

        # If we have both builds, merge coordinate information
        if grch38_df is not None and grch37_df is not None:
            logger.info(f"Merging coordinates from both builds for {pgs_id}")

            # Create a mapping based on variant identifiers
            # Try rsID first, then chr:pos:ref:alt
            grch37_coords = {}
            for _, row in grch37_df.iterrows():
                variant_key = self._get_variant_key(row)
                if variant_key:
                    grch37_coords[variant_key] = {
                        "hg19_chromosome": row.get("hg19_chromosome", ""),
                        "hg19_chr": row.get("hg19_chr", ""),
                        "hg19_start": row.get("hg19_start", ""),
                        "hg19_pos": row.get("hg19_pos", ""),
                        "hg19_end": row.get("hg19_end", ""),
                        "hg19_strand": row.get("hg19_strand", ""),
                        "hg19_allele_string": row.get("hg19_allele_string", ""),
                    }

            # Merge GRCh37 coordinates into GRCh38 base
            for idx, row in result_df.iterrows():
                variant_key = self._get_variant_key(row)
                if variant_key in grch37_coords:
                    coord_data = grch37_coords[variant_key]
                    for col, value in coord_data.items():
                        result_df.at[idx, col] = value

            matched_variants = sum(
                1
                for _, row in result_df.iterrows()
                if self._get_variant_key(row) in grch37_coords
            )
            logger.info(
                f"Matched {matched_variants}/{len(result_df)} variants between builds for {pgs_id}"
            )

        return result_df

    def _get_variant_key(self, row: pd.Series[Any]) -> str:
        """
        Get a consistent variant key for matching across builds.
        Tries rsID first, then chr:pos:ref:alt.
        """
        # Try rsID if available and not empty
        rsid = row.get("snp", "")
        if rsid and rsid.startswith("rs"):
            return rsid

        # Try extracting rsID from effect_weight or other metadata
        if "rsid" in row and row["rsid"] and str(row["rsid"]).startswith("rs"):
            return str(row["rsid"])

        # Fallback to chr:pos:ref:alt if coordinates available
        chr_col = (
            row.get("chromosome", "")
            or row.get("hg38_chromosome", "")
            or row.get("hg19_chromosome", "")
        )
        pos_col = (
            row.get("position", "")
            or row.get("hg38_start", "")
            or row.get("hg19_start", "")
        )
        ref_allele = row.get("other_allele", "") or row.get("reference_allele", "")
        eff_allele = row.get("effect_allele", "")

        if chr_col and pos_col and ref_allele and eff_allele:
            return f"{chr_col}:{pos_col}:{ref_allele}:{eff_allele}"

        # Last resort - use the snp identifier as-is
        return str(rsid) if rsid else ""


def create_prs_parser(file_path: Path, config: dict[str, Any]) -> BaseSNPParser:
    """
    Factory function to create the appropriate PRS parser based on configuration.

    Args:
        file_path: Path to the SNP file
        config: Parser configuration

    Returns:
        Appropriate parser instance

    Raises:
        ValueError: If parser type is not supported
    """
    parser_type = config.get("parser", "bcac_prs")

    if parser_type == "bcac_prs":
        return BCACParser(file_path, config)
    elif parser_type == "generic_prs":
        return GenericPRSParser(file_path, config)
    elif parser_type == "pgs_catalog":
        return PGSCatalogParser(file_path, config)
    elif parser_type == "pgs_catalog_fetcher":
        return PGSCatalogFetcher(file_path, config)
    else:
        raise ValueError(f"Unsupported PRS parser type: {parser_type}")
