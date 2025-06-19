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
from .io import create_bed_file, create_exon_bed_file
from .report_generator import ReportGenerator

logger = logging.getLogger(__name__)
console = Console()


def generate_outputs(
    df: pd.DataFrame,
    config: dict[str, Any],
    output_dir: Path,
    transcript_data: dict[str, Any] | None = None,
) -> None:
    """
    Generate all output files including Excel, CSV, BED files, and HTML report.

    Args:
        df: Final annotated DataFrame
        config: Configuration dictionary
        output_dir: Output directory path
        transcript_data: Optional transcript data for exon BED generation
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    config_manager = ConfigManager(config)

    # Generate data files in multiple formats
    _generate_data_files(df, config_manager, output_dir)

    # Generate BED files if enabled
    _generate_bed_files(df, config_manager, output_dir, transcript_data)

    # Generate HTML report if enabled
    _generate_html_report_if_enabled(df, config_manager, output_dir)


def _generate_data_files(
    df: pd.DataFrame, config_manager: ConfigManager, output_dir: Path
) -> dict[str, Path]:
    """
    Generate data files in multiple formats using the format strategy pattern.

    Args:
        df: DataFrame to save
        config_manager: Configuration manager instance
        output_dir: Output directory

    Returns:
        Dictionary mapping format names to file paths
    """
    formats = config_manager.get_output_formats()
    saver = DataFrameSaver()

    return saver.save_multiple_formats(
        df=df,
        base_path=output_dir,
        filename_base="master_panel",
        formats=formats,
        index=False,
    )


def _generate_bed_files(
    df: pd.DataFrame,
    config_manager: ConfigManager,
    output_dir: Path,
    transcript_data: dict[str, Any] | None,
) -> None:
    """
    Generate BED files if enabled in configuration.

    Args:
        df: DataFrame with gene data
        config_manager: Configuration manager instance
        output_dir: Output directory
        transcript_data: Optional transcript data for exon BED generation
    """
    # Generate germline BED file
    if config_manager.is_bed_enabled("germline") and "include" in df.columns:
        bed_path = output_dir / "germline_panel.bed"
        create_bed_file(df, bed_path, "include")

    # Generate exon BED files
    if config_manager.is_bed_enabled("exons") and "include" in df.columns:
        _generate_exon_bed_files(df, config_manager, output_dir, transcript_data)


def _generate_html_report_if_enabled(
    df: pd.DataFrame, config_manager: ConfigManager, output_dir: Path
) -> None:
    """
    Generate HTML report if enabled in configuration.

    Args:
        df: DataFrame with gene data
        config_manager: Configuration manager instance
        output_dir: Output directory
    """
    if config_manager.is_html_enabled():
        html_path = output_dir / "panel_report.html"
        _generate_html_report(df, config_manager.to_dict(), html_path)


def _generate_html_report(
    df: pd.DataFrame, config: dict[str, Any], output_path: Path
) -> None:
    """Generate HTML report using the ReportGenerator."""
    try:
        report_generator = ReportGenerator()
        report_generator.render(df, config, output_path)
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
            "[yellow]No included genes found for exon BED file generation[/yellow]"
        )
        return

    # Check if transcript data is available
    if not transcript_data:
        console.print(
            "[yellow]No transcript data available. Skipping exon BED file generation.[/yellow]"
        )
        return

    # Get configuration
    exon_padding = config_manager.get_exon_padding()
    transcript_types_config = config_manager.get_transcript_types_config()

    # Build transcript types list based on configuration
    transcript_types = _build_transcript_types_list(transcript_types_config)

    console.print(
        f"[blue]Generating exon BED files from {len(transcript_data)} genes with provided transcript data...[/blue]"
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
            f"[yellow]No {column_name} column found - skipping {transcript_type} exon BED[/yellow]"
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
            transcript_data, gene_symbol, transcript_id, transcript_type, row
        )

        if exons:
            genes_with_transcripts += 1
            all_exons.extend(exons)

    if all_exons:
        # Generate BED file for this transcript type
        bed_filename = f"exons_{transcript_type}_transcript.bed"
        bed_path = output_dir / bed_filename
        create_exon_bed_file(all_exons, bed_path, transcript_type, exon_padding)
        console.print(
            f"[green]Generated {bed_filename} with {len(all_exons)} exons from {genes_with_transcripts}/{genes_processed} genes[/green]"
        )
    else:
        console.print(
            f"[yellow]No {transcript_type} exon data available for BED generation[/yellow]"
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
        gene_data["all_transcripts"], transcript_id
    )
    if not target_transcript or "Exon" not in target_transcript:
        return []

    # Extract and process exon information
    exons = _process_transcript_exons(
        target_transcript["Exon"], gene_symbol, transcript_id, transcript_type, row
    )

    # Sort by rank to maintain exon order
    exons.sort(key=lambda x: x["rank"])
    return exons


def _find_target_transcript(
    transcripts: list[dict[str, Any]], transcript_id: str
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
