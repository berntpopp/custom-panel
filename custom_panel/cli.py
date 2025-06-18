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
    intermediate_format: str | None = typer.Option(
        None,
        "--intermediate-format",
        help="Format for intermediate files (csv, excel, parquet)",
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
    else:
        console.print(f"[red]Unknown source: {source}[/red]")
        console.print(
            "Available sources: panelapp, inhouse, acmg, manual, hpo, commercial, cosmic"
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
    config: dict[str, Any], output_manager: OutputManager | None = None
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
            from .core.io import save_panel_data

            save_panel_data(df, filepath, format_name)
            saved_files[format_name] = filepath

    # Generate BED files if enabled
    bed_config = output_config.get("bed_files", {})
    if bed_config.get("germline", False) and "include" in df.columns:
        bed_path = output_dir / "germline_panel.bed"
        create_bed_file(df, bed_path, "include")


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

    if "mane_transcript" in df.columns:
        mane_count = (~df["mane_transcript"].isna()).sum()
        table.add_row("Genes with MANE transcripts", str(mane_count))

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
