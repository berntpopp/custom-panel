"""
Command-line interface for the custom-panel tool.

This module provides the main CLI commands using Typer.
"""

import logging
import sys
import warnings
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import typer
import yaml
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

from .core.io import create_bed_file
from .core.output_manager import OutputManager
from .core.utils import normalize_count
from .engine.annotator import GeneAnnotator
from .engine.merger import PanelMerger
from .sources.a_incidental_findings import fetch_acmg_incidental_data
from .sources.b_manual_curation import fetch_manual_curation_data
from .sources.g00_inhouse_panels import fetch_inhouse_panels_data
from .sources.g01_panelapp import fetch_panelapp_data
from .sources.g02_hpo import fetch_hpo_neoplasm_data
from .sources.g03_commercial_panels import fetch_commercial_panels_data
from .sources.g04_cosmic_germline import fetch_cosmic_germline_data
from .sources.g05_clingen import fetch_clingen_data
from .sources.g06_gencc import fetch_gencc_data

# Suppress ALL deprecation warnings at module level
warnings.filterwarnings("ignore", category=DeprecationWarning)

app = typer.Typer(
    name="custom-panel",
    help="A modern Python tool for gene panel curation and aggregation from multiple genomic databases.",
    add_completion=False,
)

console = Console()
logger = logging.getLogger(__name__)


def setup_logging(log_level: str = "INFO") -> None:
    """Setup logging configuration."""
    import warnings

    # Suppress common deprecation warnings
    warnings.filterwarnings(
        "ignore", message=".*ARC4 has been moved.*", category=DeprecationWarning
    )
    warnings.filterwarnings(
        "ignore", message=".*'BaseCommand' is deprecated.*", category=DeprecationWarning
    )

    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def load_config(config_file: Optional[str] = None) -> dict[str, Any]:  # noqa: UP007
    """Load configuration from file."""
    # Always start with default config
    default_config_path = Path(__file__).parent / "config" / "default_config.yml"

    if not default_config_path.exists():
        typer.echo(
            "Default configuration file not found. Installation may be corrupted.",
            err=True,
        )
        raise typer.Exit(1)

    try:
        # Load default config first
        with open(default_config_path) as f:
            config = yaml.safe_load(f)
        typer.echo(f"Loaded default configuration from: {default_config_path}")

        # If a specific config file was provided, treat it as override
        if config_file:
            override_path = Path(config_file)
            if not override_path.exists():
                typer.echo(f"Configuration file not found: {override_path}", err=True)
                raise typer.Exit(1)

            typer.echo(f"Loading config overrides from: {override_path}")
            with open(override_path) as f:
                override_config = yaml.safe_load(f)
                if override_config:
                    config = _merge_configs(config, override_config)
        else:
            # No specific config file, check for local overrides
            local_config_path = Path("config.local.yml")
            if local_config_path.exists():
                typer.echo(f"Loading local config overrides from: {local_config_path}")
                with open(local_config_path) as f:
                    local_config = yaml.safe_load(f)
                    if local_config:
                        config = _merge_configs(config, local_config)

        return config
    except Exception as e:
        typer.echo(f"Error loading configuration: {e}", err=True)
        raise typer.Exit(1) from e


def _merge_configs(
    base_config: dict[str, Any], override_config: dict[str, Any]
) -> dict[str, Any]:
    """Recursively merge override configuration into base configuration."""
    import copy

    result = copy.deepcopy(base_config)

    for key, value in override_config.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _merge_configs(result[key], value)
        else:
            result[key] = value

    return result


@app.command()
def run(
    config_file: Optional[str] = typer.Option(  # noqa: UP007
        None, "--config-file", "-c", help="Configuration file path"
    ),
    output_dir: Optional[str] = typer.Option(  # noqa: UP007
        None, "--output-dir", "-o", help="Output directory"
    ),
    score_threshold: Optional[float] = typer.Option(  # noqa: UP007
        None,
        "--score-threshold",
        help="Override the evidence score threshold for gene inclusion",
    ),
    log_level: str = typer.Option(
        "INFO", "--log-level", help="Log level (DEBUG, INFO, WARNING, ERROR)"
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Run without generating output files"
    ),
    save_intermediate: bool = typer.Option(
        False, "--save-intermediate", help="Save intermediate files for debugging"
    ),
    intermediate_format: Optional[str] = typer.Option(
        None,
        "--intermediate-format",
        help="Format for intermediate files (csv, excel, parquet, json)",
    ),
    log_to_file: bool = typer.Option(False, "--log-to-file", help="Save logs to files"),
    structured_output: bool = typer.Option(
        True,
        "--structured-output/--flat-output",
        help="Use structured output directories",
    ),
) -> None:
    """
    Run the complete gene panel curation pipeline.
    """
    setup_logging(log_level)

    # Load configuration
    config = load_config(config_file)

    # Override configuration with command line arguments
    if output_dir:
        config.setdefault("general", {})["output_dir"] = output_dir
    if score_threshold is not None:
        config.setdefault("scoring", {}).setdefault("thresholds", {})[
            "score_threshold"
        ] = score_threshold

    # Override intermediate file and logging settings
    if save_intermediate:
        config.setdefault("output", {}).setdefault("intermediate_files", {})[
            "enabled"
        ] = True
    if intermediate_format is not None and intermediate_format in [
        "csv",
        "excel",
        "parquet",
        "json",
    ]:
        config.setdefault("output", {}).setdefault("intermediate_files", {})[
            "format"
        ] = intermediate_format
    if log_to_file:
        config.setdefault("output", {}).setdefault("file_logging", {})["enabled"] = True
    if not structured_output:
        config.setdefault("directory_structure", {})["use_structured_output"] = False

    output_dir_path = Path(config.get("general", {}).get("output_dir", "results"))

    console.print("[bold green]Starting custom-panel pipeline[/bold green]")
    console.print(f"Output directory: {output_dir_path}")

    if dry_run:
        console.print(
            "[yellow]Running in dry-run mode - no files will be generated[/yellow]"
        )

    try:
        # Initialize engines and output manager
        annotator = GeneAnnotator(config)
        merger = PanelMerger(config)
        output_manager = OutputManager(config, output_dir_path)

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            # Step 1: Fetch data from all sources
            task = progress.add_task("Fetching data from sources...", total=None)
            raw_dataframes = fetch_all_sources(config, output_manager)
            if not raw_dataframes:
                console.print("[red]No data fetched from any source. Exiting.[/red]")
                raise typer.Exit(1)
            progress.update(
                task, description=f"Fetched data from {len(raw_dataframes)} sources"
            )

            # Step 2: Centralized gene symbol standardization
            progress.update(
                task, description="Centralizing gene symbol standardization..."
            )

            # Combine all raw dataframes into one for efficient standardization
            unified_df = pd.concat(raw_dataframes, ignore_index=True)

            # Extract ALL unique symbols across all sources
            all_unique_symbols = (
                unified_df["approved_symbol"].dropna().unique().tolist()
            )
            logger.info(
                f"Total unique symbols across all sources: {len(all_unique_symbols)}"
            )

            # Make a single batch call to standardize ALL symbols
            progress.update(
                task,
                description=f"Standardizing {len(all_unique_symbols)} unique symbols...",
            )
            symbol_map = annotator.standardize_gene_symbols(all_unique_symbols)

            # Count total changes
            total_changed = sum(
                1
                for orig, info in symbol_map.items()
                if info["approved_symbol"] is not None
                and orig != info["approved_symbol"]
            )
            logger.info(
                f"Standardized {total_changed}/{len(all_unique_symbols)} gene symbols"
            )

            # Create mapping dictionaries
            approved_symbol_map = {
                k: v["approved_symbol"] for k, v in symbol_map.items()
            }
            hgnc_id_map = {
                k: v["hgnc_id"] for k, v in symbol_map.items() if v["hgnc_id"]
            }

            # Apply standardization to unified dataframe
            unified_df["approved_symbol"] = (
                unified_df["approved_symbol"]
                .map(approved_symbol_map)
                .fillna(unified_df["approved_symbol"])
            )

            # Add HGNC ID column if not exists
            if "hgnc_id" not in unified_df.columns:
                unified_df["hgnc_id"] = None

            # Apply HGNC ID mapping
            unified_df["hgnc_id"] = (
                unified_df["approved_symbol"]
                .map(lambda x: hgnc_id_map.get(x))
                .fillna(unified_df["hgnc_id"])
            )

            # Step 3: Pre-aggregate by source groups
            progress.update(task, description="Pre-aggregating source groups...")

            # Add source_group column based on configuration
            data_sources = config.get("data_sources", {})
            source_to_group = {}

            for group_name, group_config in data_sources.items():
                if isinstance(group_config, dict):
                    if group_config.get("source_group", False):
                        # This is a source group - map all its panels
                        for panel in group_config.get("panels", []):
                            if isinstance(panel, dict):
                                source_to_group[panel.get("name")] = group_name
                    else:
                        # Standalone source - it is its own group
                        source_to_group[group_name] = group_name

            # Apply source group mapping with prefix matching
            def map_source_to_group(source_name: str) -> str:
                # Try exact match first
                if source_name in source_to_group:
                    return source_to_group[source_name]

                # Try prefix matching for sources with colons (e.g., "PanelApp_UK:panel_name")
                if ":" in source_name:
                    prefix = source_name.split(":")[0]
                    if prefix in source_to_group:
                        return source_to_group[prefix]

                # Try partial matches for sources with different naming
                for key, group in source_to_group.items():
                    if key and (source_name.startswith(key) or key in source_name):
                        return group

                # If no match found, return the source name itself as a standalone group
                return source_name

            unified_df["source_group"] = unified_df["source_name"].apply(
                map_source_to_group
            )

            # Pre-aggregate each source group
            aggregated_sources = []

            for group_name, group_df in unified_df.groupby("source_group"):
                group_config = data_sources.get(group_name, {})

                # Group by gene symbol within this source group
                gene_groups = []

                for gene_symbol, gene_df in group_df.groupby("approved_symbol"):
                    # Count internal sources
                    internal_source_count = gene_df["source_name"].nunique()

                    # Calculate internal confidence score if this is a source group
                    if group_config.get("source_group", False):
                        normalization = group_config.get("normalization", {})
                        internal_confidence = normalize_count(
                            internal_source_count,
                            method=normalization.get("method", "linear"),
                            params=normalization,
                        )
                    else:
                        # Standalone sources have confidence of 1.0
                        internal_confidence = 1.0

                    # Create aggregated record for this gene in this source group
                    aggregated_record = {
                        "approved_symbol": gene_symbol,
                        "hgnc_id": gene_df["hgnc_id"].iloc[0],  # Should be same for all
                        "gene_name_reported": gene_df["gene_name_reported"].iloc[0]
                        if "gene_name_reported" in gene_df.columns
                        else gene_symbol,
                        "source_name": group_name,  # Use group name as source for aggregated data
                        "source_evidence_score": group_config.get(
                            "evidence_score", 1.0
                        ),
                        "source_details": f"{internal_source_count} sources in {str(group_name)}",
                        "source_group": group_name,
                        "internal_source_count": internal_source_count,
                        "internal_confidence_score": internal_confidence,
                        "category": group_config.get("category", "germline"),
                        "original_sources": [
                            str(s) for s in gene_df["source_name"].unique()
                        ],
                    }

                    gene_groups.append(aggregated_record)

                # Create DataFrame for this source group
                group_aggregated_df = pd.DataFrame(gene_groups)

                # Save intermediate aggregated data
                output_manager.save_standardized_data(
                    group_aggregated_df,
                    f"aggregated_{str(group_name)}",
                    {
                        g: {"approved_symbol": g, "hgnc_id": hgnc_id_map.get(g)}
                        for g in group_aggregated_df["approved_symbol"].unique()
                    },
                )

                aggregated_sources.append(group_aggregated_df)
                logger.info(
                    f"Source group '{str(group_name)}': {len(group_aggregated_df)} unique genes "
                    f"from {group_df['source_name'].nunique()} sources"
                )

            if not aggregated_sources:
                console.print(
                    "[red]No valid genes after pre-aggregation. Exiting.[/red]"
                )
                raise typer.Exit(1)

            # Step 4: Merge and score pre-aggregated sources
            progress.update(
                task, description="Merging and scoring pre-aggregated sources..."
            )
            master_df = merger.create_master_list(aggregated_sources, output_manager)
            if master_df.empty:
                console.print(
                    "[red]No genes in master list after merging. Exiting.[/red]"
                )
                raise typer.Exit(1)

            # Step 4: Annotate the final master list
            progress.update(task, description="Annotating final gene list...")
            annotated_df = annotator.annotate_genes(master_df)

            # Save annotated data
            annotation_summary = annotator.get_annotation_summary(annotated_df)
            output_manager.save_annotated_data(annotated_df, annotation_summary)

            # Step 5: Generate outputs
            if not dry_run:
                progress.update(task, description="Generating output files...")
                final_output_dir = output_manager.get_final_output_dir()
                generate_outputs(annotated_df, config, final_output_dir)

                # Save run summary
                run_summary = output_manager.get_run_summary()
                summary_path = output_manager.run_dir / "run_summary.json"
                with open(summary_path, "w") as f:
                    import json

                    json.dump(run_summary, f, indent=2)

                # Cleanup old runs if using structured output
                if config.get("directory_structure", {}).get(
                    "use_structured_output", True
                ):
                    output_manager.cleanup_old_runs()

            progress.remove_task(task)

        # Display summary
        display_summary(annotated_df, config)

        if not dry_run:
            if output_manager.use_structured:
                console.print(
                    f"[bold green]Pipeline completed successfully! Results saved to: {output_manager.run_dir}[/bold green]"
                )
                if output_manager.intermediate_enabled:
                    console.print("[blue]Intermediate files saved for debugging[/blue]")
                if output_manager.file_logging_enabled:
                    log_dir = output_manager.run_dir / output_manager.subdirs.get(
                        "logs", "logs"
                    )
                    console.print(f"[blue]Logs saved to: {log_dir}[/blue]")
            else:
                console.print(
                    f"[bold green]Pipeline completed successfully! Results saved to: {output_dir_path}[/bold green]"
                )
        else:
            console.print("[bold yellow]Dry run completed successfully![/bold yellow]")

    except Exception as e:
        console.print(f"[red]Pipeline failed: {e}[/red]")
        if log_level.upper() == "DEBUG":
            import traceback

            console.print(traceback.format_exc())
        raise typer.Exit(1) from e


@app.command()
def fetch(
    source: str = typer.Argument(
        ...,
        help="Data source to fetch (panelapp, inhouse, acmg, manual, hpo, commercial, cosmic)",
    ),
    config_file: Optional[str] = typer.Option(  # noqa: UP007
        None, "--config-file", "-c", help="Configuration file path"
    ),
    output_dir: str = typer.Option(
        "results/fetch", "--output-dir", "-o", help="Output directory"
    ),
    format: str = typer.Option(
        "parquet", "--format", "-f", help="Output format (parquet, csv, excel)"
    ),
    log_level: str = typer.Option("INFO", "--log-level", help="Log level"),
) -> None:
    """
    Fetch data from a single source.
    """
    setup_logging(log_level)
    config = load_config(config_file)

    console.print(f"[bold blue]Fetching data from: {source}[/bold blue]")

    # Fetch data based on source
    if source.lower() == "panelapp":
        df = fetch_panelapp_data(config)
    elif source.lower() == "inhouse":
        df = fetch_inhouse_panels_data(config)
    elif source.lower() == "acmg":
        df = fetch_acmg_incidental_data(config)
    elif source.lower() == "manual":
        df = fetch_manual_curation_data(config)
    elif source.lower() == "hpo":
        df = fetch_hpo_neoplasm_data(config)
    elif source.lower() == "commercial":
        df = fetch_commercial_panels_data(config)
    elif source.lower() == "cosmic":
        df = fetch_cosmic_germline_data(config)
    elif source.lower() == "clingen":
        df = fetch_clingen_data(config)
    elif source.lower() == "gencc":
        df = fetch_gencc_data(config)
    else:
        console.print(f"[red]Unknown source: {source}[/red]")
        console.print(
            "Available sources: panelapp, inhouse, acmg, manual, hpo, commercial, cosmic, clingen, gencc"
        )
        raise typer.Exit(1)

    if df.empty:
        console.print(f"[yellow]No data fetched from {source}[/yellow]")
        return

    # Save data
    output_path = Path(output_dir) / f"{source}_data.{format}"
    from .core.io import save_panel_data

    save_panel_data(df, output_path, format)

    console.print(f"[green]Saved {len(df)} records to: {output_path}[/green]")


@app.command()
def config_check(
    config_file: Optional[str] = typer.Option(  # noqa: UP007
        None, "--config-file", "-c", help="Configuration file path"
    ),
) -> None:
    """
    Validate and display configuration.
    """
    config = load_config(config_file)

    console.print("[bold blue]Configuration Validation[/bold blue]")

    # Check data sources
    data_sources = config.get("data_sources", {})

    table = Table(title="Data Sources Configuration")
    table.add_column("Source", style="cyan")
    table.add_column("Enabled", style="green")
    table.add_column("Status", style="yellow")

    for source_name, source_config in data_sources.items():
        enabled = source_config.get("enabled", True)
        status = "✓ Configured" if enabled else "✗ Disabled"

        # Add specific validation for each source
        if source_name == "inhouse_panels" and enabled:
            panels = source_config.get("panels", [])
            if not panels:
                status = "⚠ No panels configured"
            else:
                missing_files = []
                for panel in panels:
                    file_path = panel.get("file_path")
                    if file_path and not Path(file_path).exists():
                        missing_files.append(Path(file_path).name)
                if missing_files:
                    status = f"⚠ Missing files: {', '.join(missing_files)}"

        elif source_name == "manual_curation" and enabled:
            lists = source_config.get("lists", [])
            if not lists:
                status = "⚠ No lists configured"
            else:
                missing_files = []
                for list_item in lists:
                    file_path = list_item.get("file_path")
                    if file_path and not Path(file_path).exists():
                        missing_files.append(Path(file_path).name)
                if missing_files:
                    status = f"⚠ Missing files: {', '.join(missing_files)}"

        elif source_name == "hpo_neoplasm" and enabled:
            omim_file = source_config.get("omim_file_path")
            if omim_file and not Path(omim_file).exists():
                status = f"⚠ OMIM file not found: {Path(omim_file).name}"
            specific_terms = source_config.get("specific_hpo_terms", [])
            if (
                not source_config.get("use_neoplasm_search", True)
                and not specific_terms
                and not omim_file
            ):
                status = "⚠ No data sources configured"

        table.add_row(source_name, str(enabled), status)

    console.print(table)

    # Check scoring configuration
    scoring = config.get("scoring", {})
    if scoring:
        console.print("\n[bold blue]Scoring Configuration[/bold blue]")
        thresholds = scoring.get("thresholds", {})
        console.print(
            f"Score threshold: {thresholds.get('score_threshold', 'Not set')}"
        )
        console.print(f"Minimum sources: {thresholds.get('min_sources', 'Not set')}")


@app.command()
def search_panels(
    query: str = typer.Argument(..., help="Search term for panel names"),
    config_file: Optional[str] = typer.Option(  # noqa: UP007
        None, "--config-file", "-c", help="Configuration file path"
    ),
    log_level: str = typer.Option("INFO", "--log-level", help="Log level"),
) -> None:
    """
    Search for available panels in PanelApp.
    """
    setup_logging(log_level)
    config = load_config(config_file)

    from .sources.g01_panelapp import search_panelapp_panels

    console.print(f"[bold blue]Searching for panels matching: '{query}'[/bold blue]")

    panels = search_panelapp_panels(config, query)

    if not panels:
        console.print(f"[yellow]No panels found matching '{query}'[/yellow]")
        return

    table = Table(title=f"Panels matching '{query}'")
    table.add_column("ID", style="cyan")
    table.add_column("Name", style="green")
    table.add_column("Version", style="yellow")
    table.add_column("Source", style="blue")

    for panel in panels:
        table.add_row(
            str(panel.get("id", "")),
            panel.get("name", ""),
            panel.get("version", ""),
            panel.get("source", ""),
        )

    console.print(table)


def fetch_all_sources(
    config: dict[str, Any], output_manager: Optional[OutputManager] = None
) -> list[pd.DataFrame]:
    """Fetch data from all enabled sources."""
    dataframes = []

    # Define source functions - updated to match new config keys
    source_functions = {
        "PanelApp": fetch_panelapp_data,
        "Inhouse_Panels": fetch_inhouse_panels_data,
        "ACMG_Incidental_Findings": fetch_acmg_incidental_data,
        "Manual_Curation": fetch_manual_curation_data,
        "HPO_Neoplasm": fetch_hpo_neoplasm_data,
        "Commercial_Panels": fetch_commercial_panels_data,
        "COSMIC_Germline": fetch_cosmic_germline_data,
        "ClinGen": fetch_clingen_data,
        "TheGenCC": fetch_gencc_data,
    }

    data_sources = config.get("data_sources", {})

    for source_name, fetch_function in source_functions.items():
        source_config = data_sources.get(source_name, {})

        if not source_config.get("enabled", True):
            console.print(f"[yellow]Skipping disabled source: {source_name}[/yellow]")
            continue

        try:
            console.print(f"Fetching from {source_name}...")
            df = fetch_function(config)

            if not df.empty:
                # Save raw source data if output manager is available
                if output_manager:
                    output_manager.save_source_data(df, source_name)

                dataframes.append(df)
                console.print(f"[green]✓ {source_name}: {len(df)} records[/green]")
            else:
                console.print(f"[yellow]⚠ {source_name}: No data[/yellow]")

        except Exception as e:
            console.print(f"[red]✗ {source_name}: {e}[/red]")
            logging.exception(f"Error fetching from {source_name}")

    return dataframes


def _save_annotated_data_direct(df: pd.DataFrame, path: Path, format: str) -> None:
    """
    Save annotated data directly without schema validation.

    Args:
        df: DataFrame to save
        path: Output file path
        format: Output format
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    if format.lower() == "parquet":
        df.to_parquet(path, index=False, engine="pyarrow")
    elif format.lower() == "csv":
        df.to_csv(path, index=False)
    elif format.lower() == "excel":
        df.to_excel(path, index=False, engine="openpyxl")
    else:
        raise ValueError(f"Unsupported format: {format}")


def generate_outputs(
    df: pd.DataFrame, config: dict[str, Any], output_dir: Path
) -> None:
    """Generate all output files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save master panel in multiple formats
    output_config = config.get("output", {})
    formats = output_config.get("formats", ["excel", "csv", "parquet"])

    saved_files = {}
    for format_name in formats:
        if format_name in ["excel", "csv", "parquet"]:
            filepath = (
                output_dir / f"master_panel.{format_name.replace('excel', 'xlsx')}"
            )
            # Save annotated data directly without validation (has different schema)
            _save_annotated_data_direct(df, filepath, format_name)
            saved_files[format_name] = filepath

    # Generate BED files if enabled
    bed_config = output_config.get("bed_files", {})
    if bed_config.get("germline", False) and "include" in df.columns:
        bed_path = output_dir / "germline_panel.bed"
        create_bed_file(df, bed_path, "include")

    # Generate HTML report
    html_config = output_config.get("html_report", {})
    if html_config.get("enabled", True):  # Default enabled
        html_path = output_dir / "panel_report.html"
        _generate_html_report(df, config, html_path)


def _generate_html_report(df: pd.DataFrame, config: dict[str, Any], output_path: Path) -> None:
    """Generate an interactive HTML report with gene datatable and summary statistics."""
    import json
    from datetime import datetime
    
    # Calculate summary statistics
    total_genes = len(df)
    included_count = df["include"].sum() if "include" in df.columns else 0
    annotated_count = (~df["chromosome"].isna()).sum() if "chromosome" in df.columns else 0
    
    # MANE transcript counts
    mane_select_count = (~df["mane_select_transcript"].isna()).sum() if "mane_select_transcript" in df.columns else 0
    mane_clinical_count = (~df["mane_clinical_transcript"].isna()).sum() if "mane_clinical_transcript" in df.columns else 0
    
    # Score statistics
    score_stats = {}
    if "score" in df.columns:
        score_stats = {
            "mean": float(df["score"].mean()),
            "median": float(df["score"].median()),
            "min": float(df["score"].min()),
            "max": float(df["score"].max())
        }
    
    # Top scoring genes
    top_genes = []
    if "score" in df.columns and total_genes > 0:
        top_10 = df.nlargest(10, "score")
        for _, row in top_10.iterrows():
            top_genes.append({
                "gene": row["approved_symbol"],
                "score": float(row["score"])
            })
    
    # Source distribution
    source_distribution = {}
    if "source_name" in df.columns:
        source_counts = df.groupby("source_name").size().to_dict()
        source_distribution = {str(k): int(v) for k, v in source_counts.items()}
    
    # Prepare data for DataTable - include all relevant columns for toggle functionality
    all_potential_columns = [
        "approved_symbol", "hgnc_id", "gene_size", "chromosome", "gene_start", "gene_end", 
        "biotype", "gene_description", "canonical_transcript_coverage", "mane_select_coverage", 
        "mane_clinical_coverage", "score", "include", "source_count"
    ]
    
    # Default visible columns (optimized for display)
    default_visible_columns = ["approved_symbol", "gene_size", "mane_select_coverage", "score", "include", "source_count"]
    
    # Filter to only existing columns
    available_columns = [col for col in all_potential_columns if col in df.columns]
    default_visible = [col for col in default_visible_columns if col in df.columns]
    
    # Convert DataFrame to records for JSON serialization
    table_data = []
    for _, row in df.iterrows():
        record = {}
        for col in available_columns:
            value = row[col]
            # Handle NaN values and convert to JSON-serializable types
            if pd.isna(value):
                record[col] = None
            elif isinstance(value, (int, float, bool)):
                record[col] = value
            else:
                record[col] = str(value)
        
        # Add metadata for tooltips
        record["hgnc_id_tooltip"] = str(row.get("hgnc_id", "")) if not pd.isna(row.get("hgnc_id", "")) else ""
        record["source_names_tooltip"] = str(row.get("source_details", "")) if not pd.isna(row.get("source_details", "")) else ""
        
        # Calculate source count if not already present
        if "source_count" not in record or record["source_count"] is None:
            # Try to extract from source_details or use 1 as fallback
            source_details = record.get("source_names_tooltip", "")
            if "sources in" in source_details:
                try:
                    record["source_count"] = int(source_details.split(" sources in")[0])
                except:
                    record["source_count"] = 1
            else:
                record["source_count"] = 1
        
        table_data.append(record)
    
    # Prepare data for charts
    chart_data = {
        "scores": [float(row["score"]) for _, row in df.iterrows() if not pd.isna(row.get("score"))],
        "gene_sizes": [int(row["gene_size"]) for _, row in df.iterrows() if not pd.isna(row.get("gene_size"))],
        "source_counts": [record["source_count"] for record in table_data if record["source_count"] is not None],
        "transcript_sizes": []
    }
    
    # Collect transcript sizes for distribution
    for _, row in df.iterrows():
        if not pd.isna(row.get("canonical_transcript_coverage")):
            chart_data["transcript_sizes"].append(int(row["canonical_transcript_coverage"]))
        if not pd.isna(row.get("mane_select_coverage")):
            chart_data["transcript_sizes"].append(int(row["mane_select_coverage"]))
    
    # Generate HTML content
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Gene Panel Report</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        .header h1 {{
            margin: 0;
            font-size: 2.5em;
            font-weight: 300;
        }}
        .header p {{
            margin: 10px 0 0;
            font-size: 1.1em;
            opacity: 0.9;
        }}
        .content {{
            padding: 30px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: #f8f9fa;
            border: 1px solid #e9ecef;
            border-radius: 6px;
            padding: 20px;
            text-align: center;
        }}
        .stat-number {{
            font-size: 2em;
            font-weight: bold;
            color: #495057;
            margin-bottom: 5px;
        }}
        .stat-label {{
            color: #6c757d;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .section {{
            margin-bottom: 40px;
        }}
        .section h2 {{
            color: #495057;
            border-bottom: 2px solid #e9ecef;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .top-genes {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 10px;
            margin-bottom: 20px;
        }}
        .gene-card {{
            background: #fff;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            padding: 10px;
            text-align: center;
        }}
        .gene-name {{
            font-weight: bold;
            color: #495057;
        }}
        .gene-score {{
            color: #6c757d;
            font-size: 0.9em;
        }}
        #geneTable {{
            margin-top: 20px;
        }}
        .dataTables_wrapper .dataTables_length select,
        .dataTables_wrapper .dataTables_filter input {{
            border: 1px solid #ced4da;
            border-radius: 4px;
            padding: 6px 12px;
        }}
        .dataTables_wrapper .dataTables_paginate .paginate_button {{
            border: 1px solid #dee2e6;
            margin: 0 2px;
            padding: 6px 12px;
            border-radius: 4px;
        }}
        .dataTables_wrapper .dataTables_paginate .paginate_button.current {{
            background: #007bff;
            color: white !important;
            border-color: #007bff;
        }}
        table.dataTable {{
            border-collapse: collapse;
            width: 100%;
        }}
        table.dataTable th,
        table.dataTable td {{
            border: 1px solid #dee2e6;
            padding: 8px 12px;
        }}
        table.dataTable thead th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #495057;
        }}
        table.dataTable tbody tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        table.dataTable tbody tr:hover {{
            background-color: #e3f2fd;
        }}
        .score-high {{ color: #28a745; font-weight: bold; }}
        .score-medium {{ color: #ffc107; font-weight: bold; }}
        .score-low {{ color: #dc3545; font-weight: bold; }}
        .include-yes {{ color: #28a745; font-weight: bold; }}
        .include-no {{ color: #6c757d; }}
        .chart-container {{
            position: relative;
            height: 400px;
            margin-bottom: 30px;
        }}
        .charts-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 30px;
            margin-bottom: 40px;
        }}
        .chart-card {{
            background: #fff;
            border: 1px solid #dee2e6;
            border-radius: 6px;
            padding: 20px;
        }}
        .chart-title {{
            font-size: 1.2em;
            font-weight: 600;
            color: #495057;
            margin-bottom: 15px;
            text-align: center;
        }}
        .tooltip-cell {{
            cursor: help;
            position: relative;
        }}
        .tooltip-cell:hover {{
            background-color: #e3f2fd !important;
        }}
        .column-toggle {{
            display: inline-flex;
            align-items: center;
            gap: 5px;
            padding: 5px 10px;
            background: white;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            font-size: 0.9em;
        }}
        .column-toggle input[type="checkbox"] {{
            margin: 0;
        }}
        .column-toggle label {{
            margin: 0;
            cursor: pointer;
            user-select: none;
        }}
        .gene-link {{
            text-decoration: none !important;
            color: #007bff !important;
            font-weight: 500 !important;
            transition: color 0.2s ease;
        }}
        .gene-link:hover {{
            color: #0056b3 !important;
            text-decoration: underline !important;
        }}
        .gene-link:visited {{
            color: #6f42c1 !important;
        }}
    </style>
    <!-- DataTables CSS -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
    <!-- Chart.js for interactive charts -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <!-- jQuery and DataTables JS -->
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Gene Panel Report</h1>
            <p>Generated on {datetime.now().strftime("%B %d, %Y at %I:%M %p")}</p>
        </div>
        
        <div class="content">
            <!-- Summary Statistics -->
            <div class="section">
                <h2>Summary Statistics</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <div class="stat-number">{total_genes:,}</div>
                        <div class="stat-label">Total Genes</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-number">{included_count:,}</div>
                        <div class="stat-label">Panel Genes</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-number">{annotated_count:,}</div>
                        <div class="stat-label">With Coordinates</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-number">{mane_select_count:,}</div>
                        <div class="stat-label">MANE Select</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-number">{mane_clinical_count:,}</div>
                        <div class="stat-label">MANE Clinical</div>
                    </div>
                </div>
            </div>
            
            <!-- Interactive Charts -->
            <div class="section">
                <h2>Data Visualizations</h2>
                <div class="charts-grid">
                    <div class="chart-card">
                        <div class="chart-title">Score Distribution</div>
                        <div class="chart-container">
                            <canvas id="scoreChart"></canvas>
                        </div>
                    </div>
                    <div class="chart-card">
                        <div class="chart-title">Gene Size Distribution</div>
                        <div class="chart-container">
                            <canvas id="geneSizeChart"></canvas>
                        </div>
                    </div>
                    <div class="chart-card">
                        <div class="chart-title">Source Count Distribution</div>
                        <div class="chart-container">
                            <canvas id="sourceCountChart"></canvas>
                        </div>
                    </div>
                    <div class="chart-card">
                        <div class="chart-title">Transcript Size Distribution</div>
                        <div class="chart-container">
                            <canvas id="transcriptSizeChart"></canvas>
                        </div>
                    </div>
                </div>
            </div>
            
            <!-- Gene Data Table -->
            <div class="section">
                <h2>Gene Data Table</h2>
                <p>Interactive table showing all genes with sorting, filtering, and search capabilities.</p>
                
                <!-- Column Toggle Controls -->
                <div class="column-toggles" style="margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 6px;">
                    <h4 style="margin: 0 0 10px 0; color: #495057;">Toggle Columns:</h4>
                    <div id="columnToggles" style="display: flex; flex-wrap: wrap; gap: 10px;"></div>
                </div>
                
                <table id="geneTable" class="display" style="width:100%">
                    <thead>
                        <tr id="tableHeader"></tr>
                    </thead>
                </table>
            </div>
        </div>
    </div>

    <script>
        $(document).ready(function() {{
            // Gene data
            var geneData = {json.dumps(table_data)};
            var chartData = {json.dumps(chart_data)};
            var availableColumns = {json.dumps(available_columns)};
            var defaultVisible = {json.dumps(default_visible)};
            
            // Create interactive charts
            createCharts(chartData);
            
            // Initialize table with column toggles
            initializeTable(geneData, availableColumns, defaultVisible);
        }});
        
        function initializeTable(geneData, availableColumns, defaultVisible) {{
            // Column configuration
            var columnConfig = {{
                'approved_symbol': {{
                    title: 'Gene Symbol',
                    render: function(data, type, row) {{
                        if (type === 'display') {{
                            var hgncId = row.hgnc_id_tooltip || '';
                            var tooltipText = hgncId ? 'HGNC ID: ' + hgncId : 'No HGNC ID';
                            
                            // Use direct HGNC ID URL if available, otherwise fallback to search
                            var hgncUrl;
                            if (hgncId && hgncId.startsWith('HGNC:')) {{
                                hgncUrl = 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/' + hgncId;
                                tooltipText += ' (Click to view gene report)';
                            }} else {{
                                hgncUrl = 'https://www.genenames.org/tools/search/#!/?query=' + encodeURIComponent(data);
                                tooltipText += ' (Click to search in HGNC)';
                            }}
                            
                            return '<a href="' + hgncUrl + '" target="_blank" class="tooltip-cell gene-link" title="' + tooltipText + '">' + data + '</a>';
                        }}
                        return data;
                    }}
                }},
                'hgnc_id': {{ title: 'HGNC ID' }},
                'gene_size': {{
                    title: 'Gene Size (bp)',
                    render: function(data, type, row) {{
                        if (type === 'display' && data !== null) {{
                            return parseInt(data).toLocaleString();
                        }}
                        return data;
                    }}
                }},
                'chromosome': {{ title: 'Chromosome' }},
                'gene_start': {{
                    title: 'Start',
                    render: function(data, type, row) {{
                        if (type === 'display' && data !== null) {{
                            return parseInt(data).toLocaleString();
                        }}
                        return data;
                    }}
                }},
                'gene_end': {{
                    title: 'End',
                    render: function(data, type, row) {{
                        if (type === 'display' && data !== null) {{
                            return parseInt(data).toLocaleString();
                        }}
                        return data;
                    }}
                }},
                'biotype': {{ title: 'Biotype' }},
                'gene_description': {{ title: 'Description' }},
                'canonical_transcript_coverage': {{
                    title: 'Canonical Coverage (bp)',
                    render: function(data, type, row) {{
                        if (type === 'display' && data !== null) {{
                            return parseInt(data).toLocaleString();
                        }}
                        return data;
                    }}
                }},
                'mane_select_coverage': {{
                    title: 'MANE Select Coverage (bp)',
                    render: function(data, type, row) {{
                        if (type === 'display' && data !== null) {{
                            return parseInt(data).toLocaleString();
                        }}
                        return data;
                    }}
                }},
                'mane_clinical_coverage': {{
                    title: 'MANE Clinical Coverage (bp)',
                    render: function(data, type, row) {{
                        if (type === 'display' && data !== null) {{
                            return parseInt(data).toLocaleString();
                        }}
                        return data;
                    }}
                }},
                'source_count': {{
                    title: 'Sources',
                    render: function(data, type, row) {{
                        if (type === 'display') {{
                            var sourceNames = row.source_names_tooltip || 'No source details';
                            return '<span class="tooltip-cell" title="' + sourceNames + '">' + data + '</span>';
                        }}
                        return data;
                    }}
                }},
                'score': {{
                    title: 'Score',
                    render: function(data, type, row) {{
                        if (type === 'display' && data !== null && data !== undefined) {{
                            var scoreClass = data >= 3 ? 'score-high' : data >= 1.5 ? 'score-medium' : 'score-low';
                            return '<span class="' + scoreClass + '">' + parseFloat(data).toFixed(2) + '</span>';
                        }}
                        return data;
                    }}
                }},
                'include': {{
                    title: 'Included',
                    render: function(data, type, row) {{
                        if (type === 'display') {{
                            var includeClass = data ? 'include-yes' : 'include-no';
                            var includeText = data ? 'Yes' : 'No';
                            return '<span class="' + includeClass + '">' + includeText + '</span>';
                        }}
                        return data;
                    }}
                }}
            }};
            
            // Create column toggles
            var togglesContainer = $('#columnToggles');
            availableColumns.forEach(function(col) {{
                var isVisible = defaultVisible.includes(col);
                var displayName = columnConfig[col] ? columnConfig[col].title : col.replace(/_/g, ' ').replace(/\\b\\w/g, l => l.toUpperCase());
                
                var toggleHtml = '<div class="column-toggle">' +
                    '<input type="checkbox" id="toggle_' + col + '" ' + (isVisible ? 'checked' : '') + '>' +
                    '<label for="toggle_' + col + '">' + displayName + '</label>' +
                    '</div>';
                togglesContainer.append(toggleHtml);
            }});
            
            // Build initial columns
            var columns = buildColumns(availableColumns, defaultVisible, columnConfig);
            
            // Initialize DataTable
            var table = $('#geneTable').DataTable({{
                data: geneData,
                columns: columns,
                pageLength: 25,
                lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
                order: [[columns.findIndex(col => col.data === 'score'), 'desc']],
                responsive: true,
                autoWidth: false,
                scrollX: true,
                processing: true,
                language: {{
                    processing: "Loading gene data...",
                    search: "Search genes:",
                    lengthMenu: "Show _MENU_ genes per page",
                    info: "Showing _START_ to _END_ of _TOTAL_ genes",
                    infoEmpty: "Showing 0 to 0 of 0 genes",
                    infoFiltered: "(filtered from _MAX_ total genes)",
                    paginate: {{
                        first: "First",
                        last: "Last",
                        next: "Next",
                        previous: "Previous"
                    }}
                }}
            }});
            
            // Add toggle event listeners
            availableColumns.forEach(function(col) {{
                $('#toggle_' + col).on('change', function() {{
                    var column = table.column(col + ':name');
                    column.visible(this.checked);
                }});
            }});
        }}
        
        function buildColumns(availableColumns, visibleColumns, columnConfig) {{
            var columns = [];
            
            availableColumns.forEach(function(col) {{
                var config = columnConfig[col] || {{ title: col.replace(/_/g, ' ').replace(/\\b\\w/g, l => l.toUpperCase()) }};
                columns.push({{
                    data: col,
                    name: col,
                    title: config.title,
                    render: config.render || null,
                    visible: visibleColumns.includes(col)
                }});
            }});
            
            return columns;
        }}
        
        function createCharts(data) {{
            // Score Distribution with natural boundaries
            if (data.scores && data.scores.length > 0) {{
                var ctx1 = document.getElementById('scoreChart').getContext('2d');
                var scoreHistogram = createScoreHistogram(data.scores);
                new Chart(ctx1, {{
                    type: 'bar',
                    data: {{
                        labels: scoreHistogram.labels,
                        datasets: [{{
                            label: 'Gene Count',
                            data: scoreHistogram.data,
                            backgroundColor: 'rgba(102, 126, 234, 0.6)',
                            borderColor: 'rgba(102, 126, 234, 1)',
                            borderWidth: 1
                        }}]
                    }},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {{
                            legend: {{ display: false }},
                            title: {{ display: true, text: 'Distribution of Gene Scores' }}
                        }},
                        scales: {{
                            y: {{ beginAtZero: true, title: {{ display: true, text: 'Number of Genes' }} }},
                            x: {{ title: {{ display: true, text: 'Score Range' }} }}
                        }}
                    }}
                }});
            }}
            
            // Gene Size Distribution with natural boundaries
            if (data.gene_sizes && data.gene_sizes.length > 0) {{
                var ctx2 = document.getElementById('geneSizeChart').getContext('2d');
                var sizeHistogram = createSizeHistogram(data.gene_sizes);
                new Chart(ctx2, {{
                    type: 'bar',
                    data: {{
                        labels: sizeHistogram.labels,
                        datasets: [{{
                            label: 'Gene Count',
                            data: sizeHistogram.data,
                            backgroundColor: 'rgba(75, 192, 192, 0.6)',
                            borderColor: 'rgba(75, 192, 192, 1)',
                            borderWidth: 1
                        }}]
                    }},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {{
                            legend: {{ display: false }},
                            title: {{ display: true, text: 'Distribution of Gene Sizes' }}
                        }},
                        scales: {{
                            y: {{ beginAtZero: true, title: {{ display: true, text: 'Number of Genes' }} }},
                            x: {{ title: {{ display: true, text: 'Gene Size (bp)' }} }}
                        }}
                    }}
                }});
            }}
            
            // Source Count Distribution
            if (data.source_counts && data.source_counts.length > 0) {{
                var ctx3 = document.getElementById('sourceCountChart').getContext('2d');
                var sourceCountData = data.source_counts.reduce((acc, count) => {{
                    acc[count] = (acc[count] || 0) + 1;
                    return acc;
                }}, {{}});
                new Chart(ctx3, {{
                    type: 'bar',
                    data: {{
                        labels: Object.keys(sourceCountData).sort((a, b) => parseInt(a) - parseInt(b)),
                        datasets: [{{
                            label: 'Gene Count',
                            data: Object.keys(sourceCountData).sort((a, b) => parseInt(a) - parseInt(b)).map(k => sourceCountData[k]),
                            backgroundColor: 'rgba(255, 99, 132, 0.6)',
                            borderColor: 'rgba(255, 99, 132, 1)',
                            borderWidth: 1
                        }}]
                    }},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {{
                            legend: {{ display: false }},
                            title: {{ display: true, text: 'Distribution of Source Counts' }}
                        }},
                        scales: {{
                            y: {{ beginAtZero: true, title: {{ display: true, text: 'Number of Genes' }} }},
                            x: {{ title: {{ display: true, text: 'Number of Sources' }} }}
                        }}
                    }}
                }});
            }}
            
            // Transcript Size Distribution
            if (data.transcript_sizes && data.transcript_sizes.length > 0) {{
                var ctx4 = document.getElementById('transcriptSizeChart').getContext('2d');
                var transcriptHistogram = createTranscriptSizeHistogram(data.transcript_sizes);
                new Chart(ctx4, {{
                    type: 'bar',
                    data: {{
                        labels: transcriptHistogram.labels,
                        datasets: [{{
                            label: 'Transcript Count',
                            data: transcriptHistogram.data,
                            backgroundColor: 'rgba(255, 159, 64, 0.6)',
                            borderColor: 'rgba(255, 159, 64, 1)',
                            borderWidth: 1
                        }}]
                    }},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {{
                            legend: {{ display: false }},
                            title: {{ display: true, text: 'Distribution of Transcript Sizes' }}
                        }},
                        scales: {{
                            y: {{ beginAtZero: true, title: {{ display: true, text: 'Number of Transcripts' }} }},
                            x: {{ title: {{ display: true, text: 'Transcript Size (bp)' }} }}
                        }}
                    }}
                }});
            }}
        }}
        
        function createScoreHistogram(scores) {{
            // Natural boundaries for scores: 0-1, 1-2, 2-3, 3-4, 4-5, 5+
            var boundaries = [0, 1, 2, 3, 4, 5, 999];
            var labels = ['0-1', '1-2', '2-3', '3-4', '4-5', '5+'];
            var histogram = new Array(labels.length).fill(0);
            
            scores.forEach(score => {{
                for (var i = 0; i < boundaries.length - 1; i++) {{
                    if (score >= boundaries[i] && score < boundaries[i + 1]) {{
                        histogram[i]++;
                        break;
                    }}
                }}
            }});
            
            return {{ labels: labels, data: histogram }};
        }}
        
        function createSizeHistogram(sizes) {{
            // Natural boundaries for sizes: 0-10k, 10k-50k, 50k-100k, 100k-200k, 200k+
            var boundaries = [0, 10000, 50000, 100000, 200000, 500000, 999999999];
            var labels = ['0-10k', '10k-50k', '50k-100k', '100k-200k', '200k-500k', '500k+'];
            var histogram = new Array(labels.length).fill(0);
            
            sizes.forEach(size => {{
                for (var i = 0; i < boundaries.length - 1; i++) {{
                    if (size >= boundaries[i] && size < boundaries[i + 1]) {{
                        histogram[i]++;
                        break;
                    }}
                }}
            }});
            
            return {{ labels: labels, data: histogram }};
        }}
        
        function createTranscriptSizeHistogram(sizes) {{
            // More granular boundaries for transcript sizes in 2kb steps
            var boundaries = [0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000, 25000, 30000, 999999999];
            var labels = ['0-2k', '2k-4k', '4k-6k', '6k-8k', '8k-10k', '10k-12k', '12k-14k', '14k-16k', '16k-18k', '18k-20k', '20k-25k', '25k-30k', '30k+'];
            var histogram = new Array(labels.length).fill(0);
            
            sizes.forEach(size => {{
                for (var i = 0; i < boundaries.length - 1; i++) {{
                    if (size >= boundaries[i] && size < boundaries[i + 1]) {{
                        histogram[i]++;
                        break;
                    }}
                }}
            }});
            
            return {{ labels: labels, data: histogram }};
        }}
    </script>
</body>
</html>'''
    
    # Write HTML file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    console.print(f"[green]Generated interactive HTML report: {output_path}[/green]")


def display_summary(df: pd.DataFrame, config: dict[str, Any]) -> None:
    """Display pipeline summary statistics."""
    console.print("\n[bold blue]Pipeline Summary[/bold blue]")

    # Basic statistics
    table = Table(title="Gene Panel Statistics")
    table.add_column("Metric", style="cyan")
    table.add_column("Count", style="green")

    total_genes = len(df)
    table.add_row("Total genes", str(total_genes))

    if "include" in df.columns:
        included_count = df["include"].sum()
        table.add_row("Panel genes", str(included_count))

    # Annotation statistics
    if "chromosome" in df.columns:
        annotated_count = (~df["chromosome"].isna()).sum()
        table.add_row("Genes with coordinates", str(annotated_count))

    if (
        "mane_select_transcript" in df.columns
        or "mane_clinical_transcript" in df.columns
    ):
        mane_select_count = (
            (~df["mane_select_transcript"].isna()).sum()
            if "mane_select_transcript" in df.columns
            else 0
        )
        mane_clinical_count = (
            (~df["mane_clinical_transcript"].isna()).sum()
            if "mane_clinical_transcript" in df.columns
            else 0
        )
        table.add_row("Genes with MANE Select", str(mane_select_count))
        table.add_row("Genes with MANE Clinical", str(mane_clinical_count))

    console.print(table)

    # Top scoring genes
    if "score" in df.columns and total_genes > 0:
        console.print("\n[bold blue]Top 10 Scoring Genes[/bold blue]")
        top_genes = df.nlargest(10, "score")

        top_table = Table()
        top_table.add_column("Gene", style="cyan")
        top_table.add_column("Score", style="green")

        for _, row in top_genes.iterrows():
            top_table.add_row(
                row["approved_symbol"],
                f"{row['score']:.2f}",
            )

        console.print(top_table)


if __name__ == "__main__":
    app()
