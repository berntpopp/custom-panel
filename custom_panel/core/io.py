"""
Standardized I/O operations for the custom-panel package.

This module provides functions for reading and writing panel data in standardized formats,
primarily using Parquet for efficient storage and retrieval.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

# Standard DataFrame schema for panel data
STANDARD_COLUMNS = [
    "approved_symbol",  # str: HGNC approved gene symbol
    "hgnc_id",  # str: HGNC ID (e.g., "HGNC:5")
    "gene_name_reported",  # str: Original gene name from source
    "source_name",  # str: Data source identifier
    "source_evidence_score",  # float: Normalized evidence score (0.0-5.0)
    "source_details",  # str: Additional source-specific information
]

# Standard annotation columns (added by the annotation engine)
ANNOTATION_COLUMNS = [
    "gene_id",  # str: Ensembl gene ID
    "chromosome",  # str: Chromosome name
    "gene_start",  # int: Gene start position
    "gene_end",  # int: Gene end position
    "gene_strand",  # int: Gene strand (-1, 1)
    "gene_size",  # int: Gene size in base pairs
    "biotype",  # str: Gene biotype
    "gene_description",  # str: Gene description
    "canonical_transcript",  # str: Canonical transcript ID
    "gene_coverage_with_padding",  # int: Gene coverage with padding
    "mane_select_transcript",  # str: MANE Select transcript ID
    "mane_select_refseq",  # str: MANE Select RefSeq match
    "mane_clinical_transcript",  # str: MANE Plus Clinical transcript ID
    "mane_clinical_refseq",  # str: MANE Plus Clinical RefSeq match
    "canonical_transcript_coverage",  # int: Canonical transcript coverage with padding
    "mane_select_coverage",  # int: MANE Select coverage with padding
    "mane_clinical_coverage",  # int: MANE Plus Clinical coverage with padding
]


def validate_panel_dataframe(
    df: pd.DataFrame, require_annotations: bool = False
) -> bool:
    """
    Validate that a DataFrame follows the standard panel schema.

    Args:
        df: DataFrame to validate
        require_annotations: Whether to require annotation columns

    Returns:
        True if valid, False otherwise

    Raises:
        ValueError: If DataFrame is invalid
    """
    if df.empty:
        raise ValueError("DataFrame is empty")

    required_columns = STANDARD_COLUMNS.copy()
    if require_annotations:
        required_columns.extend(ANNOTATION_COLUMNS)

    missing_columns = set(required_columns) - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Validate data types
    if "source_evidence_score" in df.columns:
        if not pd.api.types.is_numeric_dtype(df["source_evidence_score"]):
            raise ValueError("source_evidence_score must be numeric")

        if (
            df["source_evidence_score"].min() < 0
            or df["source_evidence_score"].max() > 5.0
        ):
            raise ValueError("source_evidence_score must be between 0.0 and 5.0")

    return True


def create_standard_dataframe(
    genes: list[str],
    source_name: str,
    evidence_scores: list[float] | None = None,
    source_details: list[str] | None = None,
    gene_names_reported: list[str] | None = None,
) -> pd.DataFrame:
    """
    Create a standardized DataFrame from gene list and metadata.

    Args:
        genes: List of gene symbols
        source_name: Name of the data source
        evidence_scores: Evidence scores (default: 1.0 for all)
        source_details: Source-specific details
        gene_names_reported: Original gene names as reported by source

    Returns:
        Standardized DataFrame
    """
    n_genes = len(genes)

    # Set defaults
    if evidence_scores is None:
        evidence_scores = [1.0] * n_genes
    if source_details is None:
        source_details = [""] * n_genes
    if gene_names_reported is None:
        gene_names_reported = genes.copy()

    # Validate lengths
    if len(evidence_scores) != n_genes:
        raise ValueError("Length of evidence_scores must match genes")
    if len(source_details) != n_genes:
        raise ValueError("Length of source_details must match genes")
    if len(gene_names_reported) != n_genes:
        raise ValueError("Length of gene_names_reported must match genes")

    df = pd.DataFrame(
        {
            "approved_symbol": genes,  # Raw gene symbols, will be normalized by HGNC annotation
            "hgnc_id": [""] * n_genes,  # Will be filled by HGNC annotation
            "gene_name_reported": gene_names_reported,  # Original symbols from source
            "source_name": [source_name] * n_genes,
            "source_evidence_score": evidence_scores,
            "source_details": source_details,
        }
    )

    return df


def save_panel_data(
    df: pd.DataFrame, path: str | Path, format: str = "parquet"
) -> None:
    """
    Save panel data to file.

    Args:
        df: DataFrame to save
        path: Output file path
        format: Output format ("parquet", "csv", "excel", "json")

    Raises:
        ValueError: If format is not supported
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    validate_panel_dataframe(df, require_annotations=False)

    if format.lower() == "parquet":
        df.to_parquet(path, index=False, engine="pyarrow")
    elif format.lower() == "csv":
        df.to_csv(path, index=False)
    elif format.lower() == "excel":
        df.to_excel(path, index=False, engine="openpyxl")
    elif format.lower() == "json":
        df.to_json(path, orient="records", indent=2)
    else:
        raise ValueError(f"Unsupported format: {format}")

    logger.info(f"Saved {len(df)} records to {path}")


def load_panel_data(path: str | Path, format: str | None = None) -> pd.DataFrame:
    """
    Load panel data from file.

    Args:
        path: Input file path
        format: Input format (auto-detected if None)

    Returns:
        Loaded DataFrame

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If format is not supported
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # Auto-detect format if not specified
    if format is None:
        format = path.suffix.lower()[1:]  # Remove the dot

    if format == "parquet":
        df = pd.read_parquet(path, engine="pyarrow")
    elif format == "csv":
        df = pd.read_csv(path)
    elif format in ["xlsx", "xls", "excel"]:
        df = pd.read_excel(path, engine="openpyxl")
    elif format == "json":
        df = pd.read_json(path, orient="records")
    else:
        raise ValueError(f"Unsupported format: {format}")

    logger.info(f"Loaded {len(df)} records from {path}")
    return df


def save_master_panel(
    df: pd.DataFrame, output_dir: str | Path, base_name: str = "master_panel"
) -> dict[str, Path]:
    """
    Save master panel data in multiple formats.

    Args:
        df: Master panel DataFrame
        output_dir: Output directory
        base_name: Base filename (without extension)

    Returns:
        Dictionary mapping format to saved file path
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    # Save in multiple formats
    formats = [
        ("parquet", "parquet"),
        ("csv", "csv"),
        ("excel", "xlsx"),
        ("json", "json"),
    ]

    for format_name, extension in formats:
        filepath = output_dir / f"{base_name}.{extension}"
        save_panel_data(df, filepath, format_name)
        saved_files[format_name] = filepath

    return saved_files


def create_bed_file(
    df: pd.DataFrame, output_path: str | Path, filter_column: str | None = None
) -> None:
    """
    Create a BED file from annotated panel data.

    Args:
        df: Annotated DataFrame with genomic coordinates
        output_path: Output BED file path
        filter_column: Column to filter by (e.g., "include_germline")
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter data if specified
    if filter_column and filter_column in df.columns:
        df = df[df[filter_column]].copy()

    # Check for required columns
    required_cols = ["chromosome", "gene_start", "gene_end", "approved_symbol"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns for BED file: {missing_cols}")

    # Filter out rows with missing coordinates
    df = df.dropna(subset=["chromosome", "gene_start", "gene_end"])

    if df.empty:
        logger.warning(f"No genes with valid coordinates for BED file: {output_path}")
        return

    # Create BED format DataFrame
    bed_df = pd.DataFrame(
        {
            "chrom": df["chromosome"].astype(str),
            "chromStart": df["gene_start"].astype(int) - 1,  # BED is 0-based
            "chromEnd": df["gene_end"].astype(int),
            "name": df["approved_symbol"],
            "score": 1000,  # Default score
            "strand": (
                df["gene_strand"].fillna("+") if "gene_strand" in df.columns else "+"
            ),
        }
    )

    # Sort by chromosome and position with natural chromosome ordering
    # Extract numeric part of chromosome for sorting
    bed_df["chrom_num"] = (
        bed_df["chrom"]
        .str.replace("chr", "", case=False)
        .str.replace("X", "23")
        .str.replace("Y", "24")
        .str.replace("M", "25")
        .str.replace("MT", "25")  # Alternative mitochondrial notation
    )

    # Handle chromosomes that might not convert to integers (keep as high numbers)
    def safe_int_convert(x: str) -> int:
        try:
            return int(x)
        except (ValueError, TypeError):
            return 999  # Put unknown chromosomes at the end

    bed_df["chrom_num"] = bed_df["chrom_num"].apply(safe_int_convert)

    # Sort by the numeric chromosome column, then by start position
    bed_df = bed_df.sort_values(["chrom_num", "chromStart"])

    # Remove the temporary column before saving
    bed_df = bed_df.drop(columns=["chrom_num"])

    # Save BED file
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(f"Created BED file with {len(bed_df)} regions: {output_path}")


def create_exon_bed_file(
    exons_data: list[dict[str, Any]],
    output_path: str | Path,
    transcript_type: str = "canonical",
    padding: int = 0,
) -> None:
    """
    Create a BED file from exon data.

    Args:
        exons_data: List of exon dictionaries with coordinates
        output_path: Output BED file path
        transcript_type: Type of transcript (for naming)
        padding: Padding around exons in base pairs
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not exons_data:
        logger.warning(f"No exon data provided for BED file: {output_path}")
        return

    # Convert exon data to DataFrame
    exon_records = []
    for exon in exons_data:
        if all(key in exon for key in ["chromosome", "start", "end", "gene_symbol"]):
            record = {
                "chromosome": str(exon["chromosome"]),
                "start": int(exon["start"])
                - 1
                - padding,  # BED is 0-based, add padding
                "end": int(exon["end"]) + padding,
                "gene_symbol": exon["gene_symbol"],
                "exon_id": exon.get("exon_id", ""),
                "strand": exon.get("strand", "+"),
                "transcript_id": exon.get("transcript_id", ""),
                "rank": exon.get("rank", 0),
            }
            # Ensure start is not negative
            record["start"] = max(0, record["start"])
            exon_records.append(record)

    if not exon_records:
        logger.warning(f"No valid exon records for BED file: {output_path}")
        return

    exon_df = pd.DataFrame(exon_records)

    # Create BED format DataFrame
    bed_df = pd.DataFrame(
        {
            "chrom": exon_df["chromosome"],
            "chromStart": exon_df["start"],
            "chromEnd": exon_df["end"],
            "name": exon_df["gene_symbol"]
            + "_"
            + exon_df["transcript_id"]
            + "_exon"
            + exon_df["rank"].astype(str),
            "score": 1000,  # Default score
            "strand": exon_df["strand"].fillna("+"),
        }
    )

    # Sort by chromosome and position with natural chromosome ordering
    # Extract numeric part of chromosome for sorting
    bed_df["chrom_num"] = (
        bed_df["chrom"]
        .str.replace("chr", "", case=False)
        .str.replace("X", "23")
        .str.replace("Y", "24")
        .str.replace("M", "25")
        .str.replace("MT", "25")  # Alternative mitochondrial notation
    )

    # Handle chromosomes that might not convert to integers (keep as high numbers)
    def safe_int_convert(x: str) -> int:
        try:
            return int(x)
        except (ValueError, TypeError):
            return 999  # Put unknown chromosomes at the end

    bed_df["chrom_num"] = bed_df["chrom_num"].apply(safe_int_convert)

    # Sort by the numeric chromosome column, then by start position
    bed_df = bed_df.sort_values(["chrom_num", "chromStart"])

    # Remove the temporary column before saving
    bed_df = bed_df.drop(columns=["chrom_num"])

    # Save BED file
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(
        f"Created {transcript_type} exon BED file with {len(bed_df)} exons: {output_path}"
    )


def merge_panel_dataframes(dataframes: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Merge multiple panel DataFrames into a single DataFrame.

    Args:
        dataframes: List of panel DataFrames to merge

    Returns:
        Merged DataFrame
    """
    if not dataframes:
        return pd.DataFrame(columns=STANDARD_COLUMNS)

    # Validate all DataFrames
    for i, df in enumerate(dataframes):
        try:
            validate_panel_dataframe(df, require_annotations=False)
        except ValueError as e:
            logger.error(f"DataFrame {i} is invalid: {e}")
            raise

    # Concatenate all DataFrames
    merged_df = pd.concat(dataframes, ignore_index=True)

    logger.info(f"Merged {len(dataframes)} DataFrames into {len(merged_df)} records")
    return merged_df


def get_unique_genes(df: pd.DataFrame) -> list[str]:
    """
    Get unique gene symbols from a panel DataFrame.

    Args:
        df: Panel DataFrame

    Returns:
        List of unique approved gene symbols
    """
    if "approved_symbol" not in df.columns:
        raise ValueError("DataFrame must have 'approved_symbol' column")

    unique_genes = df["approved_symbol"].dropna().unique().tolist()
    logger.info(f"Found {len(unique_genes)} unique genes")
    return unique_genes


def create_genes_all_bed(
    df: pd.DataFrame, output_path: str | Path, padding: int = 0
) -> None:
    """
    Create a BED file containing all genes regardless of inclusion status.

    Args:
        df: Annotated DataFrame with genomic coordinates
        output_path: Output BED file path
        padding: Padding around genes in base pairs
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Check for required columns
    required_cols = ["chromosome", "gene_start", "gene_end", "approved_symbol"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns for BED file: {missing_cols}")

    # Filter out rows with missing coordinates
    df = df.dropna(subset=["chromosome", "gene_start", "gene_end"])

    if df.empty:
        logger.warning(f"No genes with valid coordinates for BED file: {output_path}")
        return

    # Create BED format DataFrame with padding
    bed_df = pd.DataFrame(
        {
            "chrom": df["chromosome"].astype(str),
            "chromStart": (df["gene_start"].astype(int) - 1 - padding).clip(
                lower=0
            ),  # BED is 0-based
            "chromEnd": df["gene_end"].astype(int) + padding,
            "name": df["approved_symbol"],
            "score": 1000,  # Default score
            "strand": (
                df["gene_strand"].fillna("+") if "gene_strand" in df.columns else "+"
            ),
            "element_type": "gene",
            "element_subtype": "all",
        }
    )

    # Sort by chromosome and position with natural chromosome ordering
    bed_df = _sort_bed_dataframe(bed_df)

    # Save BED file
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(f"Created genes_all BED file with {len(bed_df)} genes: {output_path}")


def create_genes_included_bed(
    df: pd.DataFrame, output_path: str | Path, padding: int = 0
) -> None:
    """
    Create a BED file containing only genes marked for inclusion.

    Args:
        df: Annotated DataFrame with genomic coordinates
        output_path: Output BED file path
        padding: Padding around genes in base pairs
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter for included genes
    if "include" not in df.columns:
        logger.warning("No 'include' column found, creating empty BED file")
        return

    df = df[df["include"]].copy()

    # Check for required columns
    required_cols = ["chromosome", "gene_start", "gene_end", "approved_symbol"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns for BED file: {missing_cols}")

    # Filter out rows with missing coordinates
    df = df.dropna(subset=["chromosome", "gene_start", "gene_end"])

    if df.empty:
        logger.warning(
            f"No included genes with valid coordinates for BED file: {output_path}"
        )
        return

    # Create BED format DataFrame with padding
    bed_df = pd.DataFrame(
        {
            "chrom": df["chromosome"].astype(str),
            "chromStart": (df["gene_start"].astype(int) - 1 - padding).clip(
                lower=0
            ),  # BED is 0-based
            "chromEnd": df["gene_end"].astype(int) + padding,
            "name": df["approved_symbol"],
            "score": 1000,  # Default score
            "strand": (
                df["gene_strand"].fillna("+") if "gene_strand" in df.columns else "+"
            ),
            "element_type": "gene",
            "element_subtype": "included",
        }
    )

    # Sort by chromosome and position with natural chromosome ordering
    bed_df = _sort_bed_dataframe(bed_df)

    # Save BED file
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(
        f"Created genes_included BED file with {len(bed_df)} genes: {output_path}"
    )


def create_snps_all_bed(
    snp_data: dict[str, pd.DataFrame], output_path: str | Path
) -> None:
    """
    Create a BED file containing all SNPs from all categories combined.

    Args:
        snp_data: Dictionary of SNP DataFrames by type
        output_path: Output BED file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not snp_data:
        logger.warning("No SNP data provided for combined BED file")
        return

    all_snps = []
    coord_columns = ["hg38_chromosome", "hg38_start", "hg38_end"]

    for snp_type, snp_df in snp_data.items():
        if snp_df.empty:
            continue

        # Check if SNP data has coordinate information
        if not all(col in snp_df.columns for col in coord_columns):
            logger.warning(f"Skipping {snp_type} SNPs - missing coordinate data")
            continue

        # Filter out SNPs without coordinates
        bed_data = snp_df.dropna(subset=coord_columns)
        if bed_data.empty:
            logger.warning(f"No {snp_type} SNPs have coordinate data")
            continue

        # Filter out rows with invalid coordinate data
        valid_coords = (
            bed_data["hg38_chromosome"].notna()
            & bed_data["hg38_start"].notna()
            & bed_data["hg38_end"].notna()
            & (bed_data["hg38_chromosome"] != "")
            & (bed_data["hg38_start"] != "")
            & (bed_data["hg38_end"] != "")
        )

        bed_data_clean = bed_data[valid_coords].copy()

        if bed_data_clean.empty:
            logger.warning(f"No {snp_type} SNPs have valid coordinate data")
            continue

        # Add SNP type information
        bed_data_clean = bed_data_clean.copy()
        bed_data_clean["snp_type"] = snp_type

        all_snps.append(bed_data_clean)

    if not all_snps:
        logger.warning("No valid SNP data found for combined BED file")
        return

    # Combine all SNP data
    combined_snps = pd.concat(all_snps, ignore_index=True)

    # Create BED format DataFrame
    bed_df = pd.DataFrame(
        {
            "chrom": "chr" + combined_snps["hg38_chromosome"].astype(str),
            "chromStart": pd.to_numeric(
                combined_snps["hg38_start"], errors="coerce"
            ).astype(int)
            - 1,  # BED is 0-based
            "chromEnd": pd.to_numeric(
                combined_snps["hg38_end"], errors="coerce"
            ).astype(int),
            "name": (
                combined_snps["snp"]
                if "snp" in combined_snps.columns
                else combined_snps.index.astype(str)
            ),
            "score": 1000,  # Default score
            "strand": "+",  # Default strand for SNPs
            "element_type": "snp",
            "element_subtype": combined_snps["snp_type"],
        }
    )

    # Sort by chromosome and position with natural chromosome ordering
    bed_df = _sort_bed_dataframe(bed_df)

    # Save BED file
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(f"Created snps_all BED file with {len(bed_df)} SNPs: {output_path}")


def create_regions_all_bed(
    regions_data: dict[str, pd.DataFrame], output_path: str | Path
) -> None:
    """
    Create a BED file containing all regions from all categories combined.

    Args:
        regions_data: Dictionary of regions DataFrames by type
        output_path: Output BED file path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not regions_data:
        logger.warning("No regions data provided for combined BED file")
        return

    all_regions = []
    coord_columns = ["chromosome", "start", "end"]

    for region_type, region_df in regions_data.items():
        if region_df.empty:
            continue

        # Check if region data has coordinate information
        if not all(col in region_df.columns for col in coord_columns):
            logger.warning(f"Skipping {region_type} regions - missing coordinate data")
            continue

        # Filter out regions without coordinates
        bed_data = region_df.dropna(subset=coord_columns)
        if bed_data.empty:
            logger.warning(f"No {region_type} regions have coordinate data")
            continue

        # Add region type information
        bed_data = bed_data.copy()
        bed_data["region_type"] = region_type

        all_regions.append(bed_data)

    if not all_regions:
        logger.warning("No valid regions data found for combined BED file")
        return

    # Combine all regions data
    combined_regions = pd.concat(all_regions, ignore_index=True)

    # Create BED format DataFrame
    bed_df = pd.DataFrame(
        {
            "chrom": combined_regions["chromosome"].astype(str),
            "chromStart": pd.to_numeric(
                combined_regions["start"], errors="coerce"
            ).astype(int)
            - 1,  # BED is 0-based
            "chromEnd": pd.to_numeric(combined_regions["end"], errors="coerce").astype(
                int
            ),
            "name": (
                combined_regions["region_name"]
                if "region_name" in combined_regions.columns
                else combined_regions.index.astype(str)
            ),
            "score": 1000,  # Default score
            "strand": "+",  # Default strand for regions
            "element_type": "region",
            "element_subtype": combined_regions["region_type"],
        }
    )

    # Sort by chromosome and position with natural chromosome ordering
    bed_df = _sort_bed_dataframe(bed_df)

    # Save BED file
    bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(
        f"Created regions_all BED file with {len(bed_df)} regions: {output_path}"
    )


def create_complete_panel_bed(
    df: pd.DataFrame,
    snp_data: dict[str, pd.DataFrame] | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
    output_path: str | Path | None = None,
    padding: int = 0,
) -> None:
    """
    Create a comprehensive BED file containing only included genes plus all SNPs and regions.

    Args:
        df: Annotated DataFrame with gene data
        snp_data: Dictionary of SNP DataFrames by type
        regions_data: Dictionary of regions DataFrames by type
        output_path: Output BED file path
        padding: Padding around genes in base pairs
    """
    if output_path is None:
        raise ValueError("output_path cannot be None")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_bed_records = []

    # Process genes - ONLY INCLUDED GENES
    if not df.empty and "include" in df.columns:
        gene_required_cols = ["chromosome", "gene_start", "gene_end", "approved_symbol"]
        gene_missing_cols = [col for col in gene_required_cols if col not in df.columns]

        if not gene_missing_cols:
            # Filter for ONLY included genes with valid coordinates
            included_genes = df[df["include"]]
            gene_df = included_genes.dropna(subset=["chromosome", "gene_start", "gene_end"])

            if not gene_df.empty:
                # Create gene BED records for included genes only
                gene_bed_df = pd.DataFrame(
                    {
                        "chrom": gene_df["chromosome"].astype(str),
                        "chromStart": (
                            gene_df["gene_start"].astype(int) - 1 - padding
                        ).clip(lower=0),
                        "chromEnd": gene_df["gene_end"].astype(int) + padding,
                        "name": gene_df["approved_symbol"],
                        "score": 1000,
                        "strand": (
                            gene_df["gene_strand"].fillna("+")
                            if "gene_strand" in gene_df.columns
                            else "+"
                        ),
                        "element_type": "gene",
                        "element_subtype": "included",
                    }
                )
                all_bed_records.append(gene_bed_df)
        else:
            logger.warning(
                f"Missing required gene columns for complete BED file: {gene_missing_cols}"
            )

    # Process SNPs
    if snp_data:
        coord_columns = ["hg38_chromosome", "hg38_start", "hg38_end"]

        for snp_type, snp_df in snp_data.items():
            if snp_df.empty:
                continue

            # Check if SNP data has coordinate information
            if not all(col in snp_df.columns for col in coord_columns):
                logger.warning(f"Skipping {snp_type} SNPs - missing coordinate data")
                continue

            # Filter out SNPs without coordinates
            bed_data = snp_df.dropna(subset=coord_columns)
            if bed_data.empty:
                continue

            # Filter out rows with invalid coordinate data
            valid_coords = (
                bed_data["hg38_chromosome"].notna()
                & bed_data["hg38_start"].notna()
                & bed_data["hg38_end"].notna()
                & (bed_data["hg38_chromosome"] != "")
                & (bed_data["hg38_start"] != "")
                & (bed_data["hg38_end"] != "")
            )

            bed_data_clean = bed_data[valid_coords].copy()

            if bed_data_clean.empty:
                continue

            # Create SNP BED records
            snp_bed_df = pd.DataFrame(
                {
                    "chrom": "chr" + bed_data_clean["hg38_chromosome"].astype(str),
                    "chromStart": pd.to_numeric(
                        bed_data_clean["hg38_start"], errors="coerce"
                    ).astype(int)
                    - 1,  # BED is 0-based
                    "chromEnd": pd.to_numeric(
                        bed_data_clean["hg38_end"], errors="coerce"
                    ).astype(int),
                    "name": (
                        bed_data_clean["snp"]
                        if "snp" in bed_data_clean.columns
                        else bed_data_clean.index.astype(str)
                    ),
                    "score": 1000,
                    "strand": "+",
                    "element_type": "snp",
                    "element_subtype": snp_type,
                }
            )
            all_bed_records.append(snp_bed_df)

    # Process regions
    if regions_data:
        coord_columns = ["chromosome", "start", "end"]

        for region_type, region_df in regions_data.items():
            if region_df.empty:
                continue

            # Check if region data has coordinate information
            if not all(col in region_df.columns for col in coord_columns):
                logger.warning(
                    f"Skipping {region_type} regions - missing coordinate data"
                )
                continue

            # Filter out regions without coordinates
            bed_data = region_df.dropna(subset=coord_columns)
            if bed_data.empty:
                continue

            # Create region BED records
            region_bed_df = pd.DataFrame(
                {
                    "chrom": bed_data["chromosome"].astype(str),
                    "chromStart": pd.to_numeric(
                        bed_data["start"], errors="coerce"
                    ).astype(int)
                    - 1,  # BED is 0-based
                    "chromEnd": pd.to_numeric(bed_data["end"], errors="coerce").astype(
                        int
                    ),
                    "name": (
                        bed_data["region_name"]
                        if "region_name" in bed_data.columns
                        else bed_data.index.astype(str)
                    ),
                    "score": 1000,
                    "strand": "+",
                    "element_type": "region",
                    "element_subtype": region_type,
                }
            )
            all_bed_records.append(region_bed_df)

    # Combine all records
    if not all_bed_records:
        logger.warning("No valid genomic elements found for complete panel BED file")
        return

    combined_bed_df = pd.concat(all_bed_records, ignore_index=True)

    # Sort by chromosome and position with natural chromosome ordering
    combined_bed_df = _sort_bed_dataframe(combined_bed_df)

    # Save BED file
    combined_bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(
        f"Created complete_panel BED file with {len(combined_bed_df)} elements: {output_path}"
    )


def create_complete_panel_exons_bed(
    df: pd.DataFrame,
    transcript_data: dict[str, Any] | None = None,
    snp_data: dict[str, pd.DataFrame] | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
    output_path: str | Path | None = None,
    padding: int = 0,
) -> None:
    """
    Create a comprehensive BED file containing exons from included genes plus all SNPs and regions.

    Args:
        df: Annotated DataFrame with gene data
        transcript_data: Dictionary containing transcript data with exon information
        snp_data: Dictionary of SNP DataFrames by type
        regions_data: Dictionary of regions DataFrames by type
        output_path: Output BED file path
        padding: Padding around exons in base pairs
    """
    if output_path is None:
        raise ValueError("output_path cannot be None")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_bed_records = []

    # Process exons from included genes
    if not df.empty and "include" in df.columns and transcript_data:
        included_genes = df[df["include"]]

        all_exons = []
        for _, row in included_genes.iterrows():
            gene_symbol = row["approved_symbol"]
            transcript_id = row.get("canonical_transcript")

            if pd.isna(transcript_id) or not transcript_id:
                continue

            # Extract exons from stored transcript data
            exons = _extract_exons_from_transcript_data(
                transcript_data, gene_symbol, transcript_id, row
            )

            if exons:
                all_exons.extend(exons)

        if all_exons:
            # Create exon BED records
            exon_bed_df = pd.DataFrame(
                {
                    "chrom": [f"chr{exon['chromosome']}" for exon in all_exons],
                    "chromStart": [max(0, exon["start"] - 1 - padding) for exon in all_exons],
                    "chromEnd": [exon["end"] + padding for exon in all_exons],
                    "name": [f"{exon['gene_symbol']}_exon{exon['rank']}" for exon in all_exons],
                    "score": 1000,
                    "strand": [exon.get("strand", "+") for exon in all_exons],
                    "element_type": "exon",
                    "element_subtype": "canonical",
                }
            )
            all_bed_records.append(exon_bed_df)

    # Add SNPs (reuse existing logic)
    if snp_data:
        coord_columns = ["hg38_chromosome", "hg38_start", "hg38_end"]

        for snp_type, snp_df in snp_data.items():
            if snp_df.empty:
                continue

            if not all(col in snp_df.columns for col in coord_columns):
                continue

            bed_data = snp_df.dropna(subset=coord_columns)
            if bed_data.empty:
                continue

            valid_coords = (
                bed_data["hg38_chromosome"].notna()
                & bed_data["hg38_start"].notna()
                & bed_data["hg38_end"].notna()
                & (bed_data["hg38_chromosome"] != "")
                & (bed_data["hg38_start"] != "")
                & (bed_data["hg38_end"] != "")
            )

            bed_data_clean = bed_data[valid_coords]
            if bed_data_clean.empty:
                continue

            snp_bed_df = pd.DataFrame(
                {
                    "chrom": "chr" + bed_data_clean["hg38_chromosome"].astype(str),
                    "chromStart": pd.to_numeric(
                        bed_data_clean["hg38_start"], errors="coerce"
                    ).astype(int) - 1,
                    "chromEnd": pd.to_numeric(
                        bed_data_clean["hg38_end"], errors="coerce"
                    ).astype(int),
                    "name": (
                        bed_data_clean["snp"]
                        if "snp" in bed_data_clean.columns
                        else bed_data_clean.index.astype(str)
                    ),
                    "score": 1000,
                    "strand": "+",
                    "element_type": "snp",
                    "element_subtype": snp_type,
                }
            )
            all_bed_records.append(snp_bed_df)

    # Add regions (reuse existing logic)
    if regions_data:
        coord_columns = ["chromosome", "start", "end"]

        for region_type, region_df in regions_data.items():
            if region_df.empty:
                continue

            if not all(col in region_df.columns for col in coord_columns):
                continue

            bed_data = region_df.dropna(subset=coord_columns)
            if bed_data.empty:
                continue

            valid_coords = (
                bed_data["chromosome"].notna()
                & bed_data["start"].notna()
                & bed_data["end"].notna()
            )

            bed_data_clean = bed_data[valid_coords]
            if bed_data_clean.empty:
                continue

            region_bed_df = pd.DataFrame(
                {
                    "chrom": "chr" + bed_data_clean["chromosome"].astype(str),
                    "chromStart": pd.to_numeric(
                        bed_data_clean["start"], errors="coerce"
                    ).astype(int) - 1,
                    "chromEnd": pd.to_numeric(
                        bed_data_clean["end"], errors="coerce"
                    ).astype(int),
                    "name": (
                        bed_data_clean["region_name"]
                        if "region_name" in bed_data_clean.columns
                        else bed_data_clean.index.astype(str)
                    ),
                    "score": 1000,
                    "strand": "+",
                    "element_type": "region",
                    "element_subtype": region_type,
                }
            )
            all_bed_records.append(region_bed_df)

    # Combine all records
    if not all_bed_records:
        logger.warning("No valid genomic elements found for complete panel exons BED file")
        return

    combined_bed_df = pd.concat(all_bed_records, ignore_index=True)
    combined_bed_df = _sort_bed_dataframe(combined_bed_df)

    # Save BED file
    combined_bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(
        f"Created complete_panel_exons BED file with {len(combined_bed_df)} elements: {output_path}"
    )


def create_complete_panel_genes_bed(
    df: pd.DataFrame,
    snp_data: dict[str, pd.DataFrame] | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
    output_path: str | Path | None = None,
    padding: int = 0,
) -> None:
    """
    Create a comprehensive BED file containing full genomic regions from included genes plus all SNPs and regions.

    Args:
        df: Annotated DataFrame with gene data
        snp_data: Dictionary of SNP DataFrames by type
        regions_data: Dictionary of regions DataFrames by type
        output_path: Output BED file path
        padding: Padding around genes in base pairs
    """
    if output_path is None:
        raise ValueError("output_path cannot be None")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    all_bed_records = []

    # Process genes - ONLY INCLUDED GENES (same as complete_panel_bed)
    if not df.empty and "include" in df.columns:
        gene_required_cols = ["chromosome", "gene_start", "gene_end", "approved_symbol"]
        gene_missing_cols = [col for col in gene_required_cols if col not in df.columns]

        if not gene_missing_cols:
            included_genes = df[df["include"]]
            gene_df = included_genes.dropna(subset=["chromosome", "gene_start", "gene_end"])

            if not gene_df.empty:
                gene_bed_df = pd.DataFrame(
                    {
                        "chrom": gene_df["chromosome"].astype(str),
                        "chromStart": (
                            gene_df["gene_start"].astype(int) - 1 - padding
                        ).clip(lower=0),
                        "chromEnd": gene_df["gene_end"].astype(int) + padding,
                        "name": gene_df["approved_symbol"],
                        "score": 1000,
                        "strand": (
                            gene_df["gene_strand"].fillna("+")
                            if "gene_strand" in gene_df.columns
                            else "+"
                        ),
                        "element_type": "gene",
                        "element_subtype": "included",
                    }
                )
                all_bed_records.append(gene_bed_df)
        else:
            logger.warning(
                f"Missing required gene columns for complete panel genes BED file: {gene_missing_cols}"
            )

    # Add SNPs (reuse existing logic)
    if snp_data:
        coord_columns = ["hg38_chromosome", "hg38_start", "hg38_end"]

        for snp_type, snp_df in snp_data.items():
            if snp_df.empty:
                continue

            if not all(col in snp_df.columns for col in coord_columns):
                continue

            bed_data = snp_df.dropna(subset=coord_columns)
            if bed_data.empty:
                continue

            valid_coords = (
                bed_data["hg38_chromosome"].notna()
                & bed_data["hg38_start"].notna()
                & bed_data["hg38_end"].notna()
                & (bed_data["hg38_chromosome"] != "")
                & (bed_data["hg38_start"] != "")
                & (bed_data["hg38_end"] != "")
            )

            bed_data_clean = bed_data[valid_coords]
            if bed_data_clean.empty:
                continue

            snp_bed_df = pd.DataFrame(
                {
                    "chrom": "chr" + bed_data_clean["hg38_chromosome"].astype(str),
                    "chromStart": pd.to_numeric(
                        bed_data_clean["hg38_start"], errors="coerce"
                    ).astype(int) - 1,
                    "chromEnd": pd.to_numeric(
                        bed_data_clean["hg38_end"], errors="coerce"
                    ).astype(int),
                    "name": (
                        bed_data_clean["snp"]
                        if "snp" in bed_data_clean.columns
                        else bed_data_clean.index.astype(str)
                    ),
                    "score": 1000,
                    "strand": "+",
                    "element_type": "snp",
                    "element_subtype": snp_type,
                }
            )
            all_bed_records.append(snp_bed_df)

    # Add regions (reuse existing logic)
    if regions_data:
        coord_columns = ["chromosome", "start", "end"]

        for region_type, region_df in regions_data.items():
            if region_df.empty:
                continue

            if not all(col in region_df.columns for col in coord_columns):
                continue

            bed_data = region_df.dropna(subset=coord_columns)
            if bed_data.empty:
                continue

            valid_coords = (
                bed_data["chromosome"].notna()
                & bed_data["start"].notna()
                & bed_data["end"].notna()
            )

            bed_data_clean = bed_data[valid_coords]
            if bed_data_clean.empty:
                continue

            region_bed_df = pd.DataFrame(
                {
                    "chrom": "chr" + bed_data_clean["chromosome"].astype(str),
                    "chromStart": pd.to_numeric(
                        bed_data_clean["start"], errors="coerce"
                    ).astype(int) - 1,
                    "chromEnd": pd.to_numeric(
                        bed_data_clean["end"], errors="coerce"
                    ).astype(int),
                    "name": (
                        bed_data_clean["region_name"]
                        if "region_name" in bed_data_clean.columns
                        else bed_data_clean.index.astype(str)
                    ),
                    "score": 1000,
                    "strand": "+",
                    "element_type": "region",
                    "element_subtype": region_type,
                }
            )
            all_bed_records.append(region_bed_df)

    # Combine all records
    if not all_bed_records:
        logger.warning("No valid genomic elements found for complete panel genes BED file")
        return

    combined_bed_df = pd.concat(all_bed_records, ignore_index=True)
    combined_bed_df = _sort_bed_dataframe(combined_bed_df)

    # Save BED file
    combined_bed_df.to_csv(output_path, sep="\t", header=False, index=False)
    logger.info(
        f"Created complete_panel_genes BED file with {len(combined_bed_df)} elements: {output_path}"
    )


def _extract_exons_from_transcript_data(
    transcript_data: dict[str, Any],
    gene_symbol: str,
    transcript_id: str,
    row: Any,
) -> list[dict[str, Any]]:
    """Extract exon data from stored transcript information."""
    gene_data = transcript_data.get(gene_symbol)
    if not gene_data or "all_transcripts" not in gene_data:
        return []

    # Find the specific transcript in the stored data
    target_transcript = None
    for transcript in gene_data["all_transcripts"]:
        if transcript.get("id") == transcript_id:
            target_transcript = transcript
            break

    if not target_transcript or "Exon" not in target_transcript:
        return []

    # Extract and process exon information
    exons = []
    for exon in target_transcript["Exon"]:
        exon_info = {
            "chromosome": exon.get("seq_region_name"),
            "start": exon.get("start"),
            "end": exon.get("end"),
            "strand": "+" if exon.get("strand", 1) == 1 else "-",
            "rank": exon.get("rank", 0),
            "gene_symbol": gene_symbol,
            "transcript_id": transcript_id,
        }
        exons.append(exon_info)

    # Sort by rank to maintain exon order
    exons.sort(key=lambda x: x["rank"])
    return exons


def _sort_bed_dataframe(bed_df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort BED DataFrame by chromosome and position with natural chromosome ordering.

    Args:
        bed_df: BED DataFrame to sort

    Returns:
        Sorted BED DataFrame
    """
    # Extract numeric part of chromosome for sorting
    bed_df["chrom_num"] = (
        bed_df["chrom"]
        .str.replace("chr", "", case=False)
        .str.replace("X", "23")
        .str.replace("Y", "24")
        .str.replace("M", "25")
        .str.replace("MT", "25")  # Alternative mitochondrial notation
    )

    # Handle chromosomes that might not convert to integers (keep as high numbers)
    def safe_int_convert(x: str) -> int:
        try:
            return int(x)
        except (ValueError, TypeError):
            return 999  # Put unknown chromosomes at the end

    bed_df["chrom_num"] = bed_df["chrom_num"].apply(safe_int_convert)

    # Sort by the numeric chromosome column, then by start position
    bed_df = bed_df.sort_values(["chrom_num", "chromStart"])

    # Remove the temporary column before returning
    bed_df = bed_df.drop(columns=["chrom_num"])

    return bed_df
