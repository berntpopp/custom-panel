"""
Improved output generator for the custom-panel tool.

This module provides functionality to generate all final output files
including Excel, CSV, Parquet, BED files, and HTML reports with better
modularization and use of DRY principles.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd
else:
    import pandas as pd

from rich.console import Console

from .config_manager import ConfigManager
from .format_strategies import DataFrameSaver
from .io import create_exon_bed_file
from .report_generator import ReportGenerator

logger = logging.getLogger(__name__)
console = Console()


def generate_outputs(
    df: pd.DataFrame,
    config: dict[str, Any],
    output_dir: Path,
    transcript_data: dict[str, Any] | None = None,
    snp_data: dict[str, pd.DataFrame] | None = None,
    deduplicated_snp_data: pd.DataFrame | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
) -> None:
    """
    Generate all output files including Excel, CSV, BED files, and HTML report.

    Args:
        df: Final annotated DataFrame
        config: Configuration dictionary
        output_dir: Output directory path
        transcript_data: Optional transcript data for exon BED generation
        snp_data: Optional SNP data dictionary by type (for individual files)
        deduplicated_snp_data: Optional deduplicated SNP data (for reports)
        regions_data: Optional regions data dictionary by type
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    config_manager = ConfigManager(config)

    # Generate data files in multiple formats (genes + SNPs + regions)
    _generate_data_files(df, config_manager, output_dir, snp_data, regions_data)

    # Generate BED files if enabled (genes + SNPs + regions)
    _generate_bed_files(
        df, config_manager, output_dir, transcript_data, snp_data, regions_data,
    )

    # Generate HTML report if enabled (genes + SNPs + regions)
    _generate_html_report_if_enabled(
        df, config_manager, output_dir, deduplicated_snp_data, regions_data,
    )


def _generate_data_files(
    df: pd.DataFrame,
    config_manager: ConfigManager,
    output_dir: Path,
    snp_data: dict[str, pd.DataFrame] | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
) -> dict[str, Path]:
    """
    Generate data files in multiple formats using the format strategy pattern.

    Args:
        df: DataFrame to save
        config_manager: Configuration manager instance
        output_dir: Output directory
        snp_data: Optional SNP data dictionary by type
        regions_data: Optional regions data dictionary by type

    Returns:
        Dictionary mapping format names to file paths
    """
    formats = config_manager.get_output_formats()
    saver = DataFrameSaver()

    # Generate gene panel files (now including regions in Excel)
    gene_files = saver.save_multiple_formats(
        df=df,
        base_path=output_dir,
        filename_base="master_panel",
        formats=formats,
        index=False,
        snp_data=snp_data,  # Pass SNP data for multi-sheet Excel
        regions_data=regions_data,  # Pass regions data for multi-sheet Excel
    )

    # Generate individual SNP files if SNP data exists
    if snp_data:
        _generate_snp_data_files(snp_data, saver, formats, output_dir)

    # Generate individual regions files if regions data exists
    if regions_data:
        _generate_regions_data_files(regions_data, saver, formats, output_dir)

    return gene_files


def _generate_snp_data_files(
    snp_data: dict[str, pd.DataFrame],
    saver: DataFrameSaver,
    formats: list[str],
    output_dir: Path,
) -> None:
    """
    Generate individual SNP data files for each SNP type.

    Args:
        snp_data: Dictionary of SNP DataFrames by type
        saver: DataFrameSaver instance
        formats: List of output formats
        output_dir: Output directory
    """
    # Combine all SNPs into a master SNP table
    all_snps = []
    for _, snp_df in snp_data.items():
        if not snp_df.empty:
            # Note: snp_type info already preserved in category column
            snp_df_copy = snp_df.copy()
            all_snps.append(snp_df_copy)

    if all_snps:
        master_snp_df = pd.concat(all_snps, ignore_index=True)

        # Clean coordinate columns for parquet compatibility
        master_snp_df = _clean_coordinate_columns(master_snp_df)

        # Generate master SNP files
        console.print(
            f"[blue]Generating master SNP files with {len(master_snp_df)} SNPs...[/blue]",
        )
        saver.save_multiple_formats(
            df=master_snp_df,
            base_path=output_dir,
            filename_base="master_snps",
            formats=formats,
            index=False,
        )

        # Generate individual SNP type files
        for snp_type, snp_df in snp_data.items():
            if not snp_df.empty:
                # Clean coordinate columns for parquet compatibility
                snp_df_clean = _clean_coordinate_columns(snp_df)
                console.print(
                    f"[blue]Generating {snp_type} SNP files with {len(snp_df_clean)} SNPs...[/blue]",
                )
                saver.save_multiple_formats(
                    df=snp_df_clean,
                    base_path=output_dir,
                    filename_base=f"snps_{snp_type}",
                    formats=formats,
                    index=False,
                )


def _generate_regions_data_files(
    regions_data: dict[str, pd.DataFrame],
    saver: DataFrameSaver,
    formats: list[str],
    output_dir: Path,
) -> None:
    """
    Generate individual regions data files for each region type.

    Args:
        regions_data: Dictionary of regions DataFrames by type
        saver: DataFrameSaver instance
        formats: List of output formats
        output_dir: Output directory
    """
    # Combine all regions into a master regions table
    all_regions = []
    for region_type, region_df in regions_data.items():
        if not region_df.empty:
            # Add region type column for identification
            region_df_copy = region_df.copy()
            region_df_copy["region_type"] = region_type
            all_regions.append(region_df_copy)

    if all_regions:
        master_regions_df = pd.concat(all_regions, ignore_index=True)
        console.print(
            f"[blue]Generating master regions files with {len(master_regions_df)} regions...[/blue]",
        )
        saver.save_multiple_formats(
            df=master_regions_df,
            base_path=output_dir,
            filename_base="master_regions",
            formats=formats,
            index=False,
        )

        # Generate individual region type files
        for region_type, region_df in regions_data.items():
            if not region_df.empty:
                console.print(
                    f"[blue]Generating {region_type} regions files with {len(region_df)} regions...[/blue]",
                )
                saver.save_multiple_formats(
                    df=region_df,
                    base_path=output_dir,
                    filename_base=f"regions_{region_type}",
                    formats=formats,
                    index=False,
                )


def _deduplicate_snps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate SNPs based on VCF ID while preserving source information.

    When the same genomic variant (VCF ID) appears in multiple sources,
    this function consolidates them into a single row with merged metadata.

    Args:
        df: DataFrame with potentially duplicate SNPs

    Returns:
        DataFrame with deduplicated SNPs and consolidated source information
    """
    if df.empty or "snp" not in df.columns:
        return df

    # Count initial entries
    initial_count = len(df)

    # Group by SNP (VCF ID) and aggregate
    aggregation_rules = {
        # Core identification - take first non-null value
        "rsid": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        # Source information - merge all sources
        "source": lambda x: "; ".join(x.dropna().unique()),
        "category": lambda x: "; ".join(x.dropna().unique()),
        # Coordinate information - take first valid coordinate
        "hg38_chromosome": lambda x: x.dropna().iloc[0]
        if not x.dropna().empty
        else None,
        "hg38_start": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "hg38_end": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "hg38_strand": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        # Allele information - take first valid allele
        "ref_allele": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "alt_allele": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "effect_allele": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "other_allele": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        # Metadata - take first valid value or merge as appropriate
        "effect_weight": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "pgs_id": lambda x: "; ".join(x.dropna().unique()),
        "pgs_name": lambda x: "; ".join(x.dropna().unique()),
        "trait": lambda x: "; ".join(x.dropna().unique()),
        "pmid": lambda x: "; ".join(x.dropna().unique()),
        "gene": lambda x: "; ".join(x.dropna().unique()),
        "gene_id": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "variant_id": lambda x: x.dropna().iloc[0] if not x.dropna().empty else None,
        "clinical_significance": lambda x: "; ".join(x.dropna().unique()),
        "review_status": lambda x: "; ".join(x.dropna().unique()),
        "consequence": lambda x: "; ".join(x.dropna().unique()),
        "distance_to_exon": lambda x: x.dropna().iloc[0]
        if not x.dropna().empty
        else None,
        "panel_name": lambda x: "; ".join(x.dropna().unique()),
        "panel_description": lambda x: "; ".join(x.dropna().unique()),
    }

    # Only apply aggregation rules to columns that exist in the DataFrame
    final_agg_rules = {
        col: rule for col, rule in aggregation_rules.items() if col in df.columns
    }

    # Group by SNP and aggregate
    deduplicated_df = df.groupby("snp", as_index=False).agg(final_agg_rules)

    # Add source count information
    source_counts = df.groupby("snp")["source"].apply(
        lambda x: len(x.dropna().unique()),
    )
    deduplicated_df["source_count"] = deduplicated_df["snp"].map(source_counts)

    # Log deduplication results
    final_count = len(deduplicated_df)
    duplicates_removed = initial_count - final_count

    if duplicates_removed > 0:
        logger.info(
            f"SNP deduplication: {duplicates_removed} duplicate entries removed "
            f"({initial_count} -> {final_count} unique SNPs)",
        )
        console.print(
            f"[yellow]Deduplicated {duplicates_removed} SNP entries - "
            f"{final_count} unique variants from {initial_count} total entries[/yellow]",
        )
    else:
        logger.info(f"No duplicate SNPs found - {final_count} unique variants")

    return deduplicated_df


def _clean_coordinate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean coordinate columns to ensure they're properly formatted for parquet.

    Converts empty strings to NaN and ensures coordinate columns are numeric.

    Args:
        df: DataFrame with potential coordinate columns

    Returns:
        DataFrame with cleaned coordinate columns
    """
    df_cleaned = df.copy()

    # Define coordinate columns that should be numeric
    coordinate_cols = [
        "hg38_start",
        "hg38_end",
        "hg38_pos",
        # Only hg38 coordinates are supported
        "start",
        "end",
        "pos",
        "position",
        "chromosome",
        "chr",
        "hm_pos",
    ]

    for col in coordinate_cols:
        if col in df_cleaned.columns:
            # Convert empty strings to NaN, then to numeric
            df_cleaned[col] = df_cleaned[col].replace("", pd.NA)
            if col in [
                "hg38_start",
                "hg38_end",
                # Only hg38 coordinates are supported
                "hm_pos",
                "start",
                "end",
                "pos",
                "position",
            ]:
                # These should be integers (positions)
                df_cleaned[col] = pd.to_numeric(
                    df_cleaned[col], errors="coerce",
                ).astype("Int64")
            elif col in ["chromosome", "chr", "hg38_chromosome"]:
                # These should be strings but handle NaN properly
                df_cleaned[col] = df_cleaned[col].astype("string")

    return df_cleaned


def _generate_bed_files(
    df: pd.DataFrame,
    config_manager: ConfigManager,
    output_dir: Path,
    transcript_data: dict[str, Any] | None,
    snp_data: dict[str, pd.DataFrame] | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
) -> None:
    """
    Generate BED files if enabled in configuration.

    Args:
        df: DataFrame with gene data
        config_manager: Configuration manager instance
        output_dir: Output directory
        transcript_data: Optional transcript data for exon BED generation
        snp_data: Optional SNP data dictionary by type
        regions_data: Optional regions data dictionary by type
    """
    # Import the new BED file functions
    from custom_panel.core.io import (
        create_complete_panel_bed,
        create_complete_panel_exons_bed,
        create_complete_panel_genes_bed,
        create_genes_all_bed,
        create_genes_included_bed,
        create_regions_all_bed,
        create_snps_all_bed,
    )

    # Get BED padding from config
    padding = config_manager.get_bed_padding()

    # Generate complete panel BED file (included genes + SNPs + regions)
    if config_manager.is_bed_enabled("complete_panel"):
        bed_path = output_dir / "complete_panel.bed"
        create_complete_panel_bed(df, snp_data, regions_data, bed_path, padding)

    # Generate complete panel exons BED file (exons from included genes + SNPs + regions)
    if config_manager.is_bed_enabled("complete_panel_exons"):
        bed_path = output_dir / "complete_panel_exons.bed"
        create_complete_panel_exons_bed(
            df, transcript_data, snp_data, regions_data, bed_path, padding,
        )

    # Generate complete panel genes BED file (full genomic regions from included genes + SNPs + regions)
    if config_manager.is_bed_enabled("complete_panel_genes"):
        bed_path = output_dir / "complete_panel_genes.bed"
        create_complete_panel_genes_bed(df, snp_data, regions_data, bed_path, padding)

    # Generate genes all BED file (all genes regardless of inclusion)
    if config_manager.is_bed_enabled("genes_all"):
        bed_path = output_dir / "genes_all.bed"
        create_genes_all_bed(df, bed_path, padding)

    # Generate genes included BED file (only included genes)
    if config_manager.is_bed_enabled("genes_included") and "include" in df.columns:
        bed_path = output_dir / "genes_included.bed"
        create_genes_included_bed(df, bed_path, padding)

    # Generate all SNPs combined BED file
    if config_manager.is_bed_enabled("snps_all") and snp_data:
        bed_path = output_dir / "snps_all.bed"
        create_snps_all_bed(snp_data, bed_path)

    # Generate all regions combined BED file
    if config_manager.is_bed_enabled("regions_all") and regions_data:
        bed_path = output_dir / "regions_all.bed"
        create_regions_all_bed(regions_data, bed_path)

    # Generate individual category BED files if enabled
    if config_manager.is_bed_enabled("individual_categories"):
        # Generate exon BED files
        if config_manager.is_bed_enabled("exons") and "include" in df.columns:
            _generate_exon_bed_files(df, config_manager, output_dir, transcript_data)

        # Generate individual SNP BED files if SNP data exists
        if snp_data:
            _generate_snp_bed_files(snp_data, output_dir)

        # Generate individual regions BED files if regions data exists
        if regions_data:
            _generate_regions_bed_files(regions_data, output_dir)

    # Legacy support: Generate genes included BED file (previously "germline_panel.bed")
    if config_manager.is_bed_enabled("germline") and "include" in df.columns:
        # Only generate if genes_included is not already enabled (avoid duplicates)
        if not config_manager.is_bed_enabled("genes_included"):
            bed_path = output_dir / "genes_included.bed"
            create_genes_included_bed(df, bed_path, padding)


def _generate_snp_bed_files(
    snp_data: dict[str, pd.DataFrame],
    output_dir: Path,
) -> None:
    """
    Generate BED files for SNP data.

    Args:
        snp_data: Dictionary of SNP DataFrames by type
        output_dir: Output directory
    """
    for snp_type, snp_df in snp_data.items():
        if snp_df.empty:
            continue

        # Check if SNP data has coordinate information
        coord_columns = ["hg38_chromosome", "hg38_start", "hg38_end"]
        if not all(col in snp_df.columns for col in coord_columns):
            console.print(
                f"[yellow]Skipping BED file for {snp_type} SNPs - missing coordinate data[/yellow]",
            )
            continue

        # Filter out SNPs without coordinates
        bed_data = snp_df.dropna(subset=coord_columns)
        if bed_data.empty:
            console.print(
                f"[yellow]No {snp_type} SNPs have coordinate data for BED file[/yellow]",
            )
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
            console.print(
                f"[yellow]No {snp_type} SNPs have valid coordinate data for BED file[/yellow]",
            )
            continue

        # Create BED format DataFrame
        bed_df = pd.DataFrame(
            {
                "chrom": "chr" + bed_data_clean["hg38_chromosome"].astype(str),
                "chromStart": pd.to_numeric(
                    bed_data_clean["hg38_start"], errors="coerce",
                ).astype(int)
                - 1,  # BED is 0-based
                "chromEnd": pd.to_numeric(
                    bed_data_clean["hg38_end"], errors="coerce",
                ).astype(int),
                "name": (
                    bed_data_clean["snp"]
                    if "snp" in bed_data_clean.columns
                    else bed_data_clean.index
                ),
                "score": 1000,  # Default score
                "strand": ".",  # Unknown strand for SNPs
            },
        )

        # Sort by chromosome and position
        bed_df = bed_df.sort_values(["chrom", "chromStart"])

        # Save BED file with improved naming
        # Map SNP types to clearer names
        snp_type_mapping = {
            "deep_intronic_clinvar": "snps_clinvar_deep_intronic",
            "identity": "snps_identity",
            "ethnicity": "snps_ethnicity",
            "prs": "snps_prs",
            "manual_snps": "snps_manual",
        }

        bed_filename = snp_type_mapping.get(snp_type, f"snps_{snp_type}")
        bed_path = output_dir / f"{bed_filename}.bed"
        bed_df.to_csv(bed_path, sep="\t", header=False, index=False)
        console.print(
            f"[green]Generated {bed_path.name} with {len(bed_df)} SNPs[/green]",
        )


def _generate_regions_bed_files(
    regions_data: dict[str, pd.DataFrame],
    output_dir: Path,
) -> None:
    """
    Generate BED files for regions data.

    Args:
        regions_data: Dictionary of regions DataFrames by type
        output_dir: Output directory
    """
    for region_type, region_df in regions_data.items():
        if region_df.empty:
            continue

        # Check if regions data has coordinate information
        coord_columns = ["chromosome", "start", "end"]
        if not all(col in region_df.columns for col in coord_columns):
            console.print(
                f"[yellow]Skipping BED file for {region_type} regions - missing coordinate data[/yellow]",
            )
            continue

        # Filter out regions without coordinates
        bed_data = region_df.dropna(subset=coord_columns)
        if bed_data.empty:
            console.print(
                f"[yellow]No {region_type} regions have coordinate data for BED file[/yellow]",
            )
            continue

        # Filter out rows with invalid coordinate data
        valid_coords = (
            bed_data["chromosome"].notna()
            & bed_data["start"].notna()
            & bed_data["end"].notna()
        )
        bed_data_clean = bed_data[valid_coords].copy()

        if bed_data_clean.empty:
            console.print(
                f"[yellow]No valid {region_type} regions for BED file[/yellow]",
            )
            continue

        # Create BED DataFrame
        bed_df = pd.DataFrame(
            {
                "chrom": "chr" + bed_data_clean["chromosome"].astype(str),
                "chromStart": bed_data_clean["start"].astype(int),
                "chromEnd": bed_data_clean["end"].astype(int),
                "name": (
                    bed_data_clean["region_name"]
                    if "region_name" in bed_data_clean.columns
                    else bed_data_clean.index
                ),
                "score": 1000,  # Default score
                "strand": ".",  # Unknown strand for regions
            },
        )

        # Sort by chromosome and position
        bed_df = bed_df.sort_values(["chrom", "chromStart"])

        # Save BED file
        bed_path = output_dir / f"regions_{region_type}.bed"
        bed_df.to_csv(bed_path, sep="\t", header=False, index=False)
        console.print(
            f"[green]Generated {bed_path.name} with {len(bed_df)} regions[/green]",
        )


def _generate_html_report_if_enabled(
    df: pd.DataFrame,
    config_manager: ConfigManager,
    output_dir: Path,
    deduplicated_snp_data: pd.DataFrame | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
) -> None:
    """
    Generate HTML report if enabled in configuration.

    Args:
        df: DataFrame with gene data
        config_manager: Configuration manager instance
        output_dir: Output directory
        deduplicated_snp_data: Optional deduplicated SNP data
        regions_data: Optional regions data dictionary by type
    """
    if config_manager.is_html_enabled():
        html_path = output_dir / "panel_report.html"
        _generate_html_report(
            df, config_manager.to_dict(), html_path, deduplicated_snp_data, regions_data,
        )


def _generate_html_report(
    df: pd.DataFrame,
    config: dict[str, Any],
    output_path: Path,
    deduplicated_snp_data: pd.DataFrame | None = None,
    regions_data: dict[str, pd.DataFrame] | None = None,
) -> None:
    """Generate HTML report using the ReportGenerator."""
    try:
        report_generator = ReportGenerator()
        report_generator.render(df, config, output_path, deduplicated_snp_data, regions_data)
    except Exception as e:
        logger.error(f"Failed to generate HTML report: {e}")
        console.print(f"[red]Failed to generate HTML report: {e}[/red]")


def _generate_exon_bed_files(
    df: pd.DataFrame,
    config_manager: ConfigManager,
    output_dir: Path,
    transcript_data: dict[str, Any] | None,
) -> None:
    """Generate exon BED files using stored transcript data from annotation."""
    # Get only included genes
    included_df = df[df["include"]].copy()
    if included_df.empty:
        console.print(
            "[yellow]No included genes found for exon BED file generation[/yellow]",
        )
        return

    # Check if transcript data is available
    if not transcript_data:
        console.print(
            "[yellow]No transcript data available. Skipping exon BED file generation.[/yellow]",
        )
        return

    # Get configuration
    exon_padding = config_manager.get_exon_padding()
    transcript_types_config = config_manager.get_transcript_types_config()

    # Build transcript types list based on configuration
    transcript_types = _build_transcript_types_list(transcript_types_config)

    console.print(
        f"[blue]Generating exon BED files from {len(transcript_data)} genes with provided transcript data...[/blue]",
    )

    # Generate BED file for each transcript type
    for transcript_type, column_name in transcript_types:
        _generate_single_transcript_bed_file(
            included_df,
            transcript_data,
            transcript_type,
            column_name,
            output_dir,
            exon_padding,
        )


def _build_transcript_types_list(
    transcript_config: dict[str, bool],
) -> list[tuple[str, str]]:
    """
    Build list of transcript types to process based on configuration.

    Args:
        transcript_config: Dictionary mapping transcript types to enabled status

    Returns:
        List of (transcript_type, column_name) tuples
    """
    transcript_types = []

    if transcript_config.get("canonical", True):
        transcript_types.append(("canonical", "canonical_transcript"))
    if transcript_config.get("mane_select", True):
        transcript_types.append(("mane_select", "mane_select_transcript"))
    if transcript_config.get("mane_clinical", False):
        transcript_types.append(("mane_clinical", "mane_clinical_transcript"))

    return transcript_types


def _generate_single_transcript_bed_file(
    included_df: pd.DataFrame,
    transcript_data: dict[str, Any],
    transcript_type: str,
    column_name: str,
    output_dir: Path,
    exon_padding: int,
) -> None:
    """
    Generate BED file for a single transcript type.

    Args:
        included_df: DataFrame with included genes
        transcript_data: Dictionary containing transcript data
        transcript_type: Type of transcript (canonical, mane_select, etc.)
        column_name: Column name in DataFrame
        output_dir: Output directory
        exon_padding: Padding to add to exons
    """
    if column_name not in included_df.columns:
        console.print(
            f"[yellow]No {column_name} column found - skipping {transcript_type} exon BED[/yellow]",
        )
        return

    all_exons = []
    genes_with_transcripts = 0
    genes_processed = 0

    for _, row in included_df.iterrows():
        gene_symbol = row["approved_symbol"]
        transcript_id = row.get(column_name)

        if pd.isna(transcript_id) or not transcript_id:
            continue

        genes_processed += 1

        # Extract exons from stored transcript data
        exons = _extract_exons_from_stored_data(
            transcript_data, gene_symbol, transcript_id, transcript_type, row,
        )

        if exons:
            genes_with_transcripts += 1
            all_exons.extend(exons)

    if all_exons:
        # Generate BED file for this transcript type with improved naming
        bed_filename = f"genes_exons_{transcript_type}.bed"
        bed_path = output_dir / bed_filename
        create_exon_bed_file(all_exons, bed_path, transcript_type, exon_padding)
        console.print(
            f"[green]Generated {bed_filename} with {len(all_exons)} exons from {genes_with_transcripts}/{genes_processed} genes[/green]",
        )
    else:
        console.print(
            f"[yellow]No {transcript_type} exon data available for BED generation[/yellow]",
        )


def _extract_exons_from_stored_data(
    transcript_data: dict[str, Any],
    gene_symbol: str,
    transcript_id: str,
    transcript_type: str,
    row: pd.Series[Any],
) -> list[dict[str, Any]]:
    """Extract exon data from stored transcript information."""
    # Get stored transcript data for this gene
    gene_data = transcript_data.get(gene_symbol)
    if not gene_data or "all_transcripts" not in gene_data:
        return []

    # Find the specific transcript in the stored data
    target_transcript = _find_target_transcript(
        gene_data["all_transcripts"], transcript_id,
    )
    if not target_transcript or "Exon" not in target_transcript:
        return []

    # Extract and process exon information
    exons = _process_transcript_exons(
        target_transcript["Exon"], gene_symbol, transcript_id, transcript_type, row,
    )

    # Sort by rank to maintain exon order
    exons.sort(key=lambda x: x["rank"])
    return exons


def _find_target_transcript(
    transcripts: list[dict[str, Any]], transcript_id: str,
) -> dict[str, Any] | None:
    """Find target transcript by ID in the transcripts list."""
    for transcript in transcripts:
        if transcript.get("id") == transcript_id:
            return transcript
    return None


def _process_transcript_exons(
    exons_data: list[dict[str, Any]],
    gene_symbol: str,
    transcript_id: str,
    transcript_type: str,
    row: pd.Series[Any],
) -> list[dict[str, Any]]:
    """Process exon data into standardized format."""
    exons = []

    for exon in exons_data:
        exon_info = {
            "exon_id": exon.get("id"),
            "chromosome": exon.get("seq_region_name"),
            "start": exon.get("start"),
            "end": exon.get("end"),
            "strand": exon.get("strand"),
            "rank": exon.get("rank", 0),
            "gene_symbol": gene_symbol,
            "gene_id": row.get("gene_id", ""),
            "transcript_id": transcript_id,
            "transcript_type": transcript_type,
        }
        exons.append(exon_info)

    return exons
