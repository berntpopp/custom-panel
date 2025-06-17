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
from .engine.annotator import GeneAnnotator
from .engine.merger import PanelMerger
from .sources.a_incidental_findings import fetch_acmg_incidental_data
from .sources.b_manual_curation import fetch_manual_curation_data
from .sources.g00_inhouse_panels import fetch_inhouse_panels_data
from .sources.g01_panelapp import fetch_panelapp_data
from .sources.g02_hpo import fetch_hpo_neoplasm_data
from .sources.g03_commercial_panels import fetch_commercial_panels_data
from .sources.s01_cosmic import fetch_cosmic_data
from .sources.s02_somatic_commercial import fetch_somatic_commercial_data

# Suppress ALL deprecation warnings at module level
warnings.filterwarnings("ignore", category=DeprecationWarning)

app = typer.Typer(
    name="custom-panel",
    help="A modern Python tool for gene panel curation and aggregation from multiple genomic databases.",
    add_completion=False,
)

console = Console()


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
    germline_threshold: Optional[float] = typer.Option(  # noqa: UP007
        None, "--germline-threshold", help="Override germline score threshold"
    ),
    somatic_threshold: Optional[float] = typer.Option(  # noqa: UP007
        None, "--somatic-threshold", help="Override somatic score threshold"
    ),
    log_level: str = typer.Option(
        "INFO", "--log-level", help="Log level (DEBUG, INFO, WARNING, ERROR)"
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Run without generating output files"
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
    if germline_threshold is not None:
        config.setdefault("scoring", {}).setdefault("thresholds", {})[
            "germline_threshold"
        ] = germline_threshold
    if somatic_threshold is not None:
        config.setdefault("scoring", {}).setdefault("thresholds", {})[
            "somatic_threshold"
        ] = somatic_threshold

    output_dir_path = Path(config.get("general", {}).get("output_dir", "results"))

    console.print("[bold green]Starting custom-panel pipeline[/bold green]")
    console.print(f"Output directory: {output_dir_path}")

    if dry_run:
        console.print(
            "[yellow]Running in dry-run mode - no files will be generated[/yellow]"
        )

    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            # Step 1: Fetch data from all sources
            task = progress.add_task("Fetching data from sources...", total=None)
            dataframes = fetch_all_sources(config)

            if not dataframes:
                console.print("[red]No data fetched from any source. Exiting.[/red]")
                raise typer.Exit(1)

            progress.update(
                task, description=f"Fetched data from {len(dataframes)} sources"
            )

            # Step 2: Merge and score genes
            progress.update(task, description="Merging and scoring genes...")
            merger = PanelMerger(config)
            master_df = merger.create_master_list(dataframes)

            if master_df.empty:
                console.print("[red]No genes in master list. Exiting.[/red]")
                raise typer.Exit(1)

            # Step 3: Annotate genes
            progress.update(task, description="Annotating genes...")
            annotator = GeneAnnotator(config)
            annotated_df = annotator.annotate_genes(master_df)

            # Step 4: Generate outputs
            if not dry_run:
                progress.update(task, description="Generating output files...")
                generate_outputs(annotated_df, config, output_dir_path)

            progress.remove_task(task)

        # Display summary
        display_summary(annotated_df, config)

        if not dry_run:
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
        help="Data source to fetch (panelapp, inhouse, acmg, manual, hpo, commercial, cosmic, somatic_commercial)",
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
        df = fetch_cosmic_data(config)
    elif source.lower() == "somatic_commercial":
        df = fetch_somatic_commercial_data(config)
    else:
        console.print(f"[red]Unknown source: {source}[/red]")
        console.print(
            "Available sources: panelapp, inhouse, acmg, manual, hpo, commercial, cosmic, somatic_commercial"
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
            f"Germline threshold: {thresholds.get('germline_threshold', 'Not set')}"
        )
        console.print(
            f"Somatic threshold: {thresholds.get('somatic_threshold', 'Not set')}"
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


def fetch_all_sources(config: dict[str, Any]) -> list[pd.DataFrame]:
    """Fetch data from all enabled sources."""
    dataframes = []

    # Define source functions
    source_functions = {
        "panelapp": fetch_panelapp_data,
        "inhouse_panels": fetch_inhouse_panels_data,
        "acmg_incidental": fetch_acmg_incidental_data,
        "manual_curation": fetch_manual_curation_data,
        "hpo_neoplasm": fetch_hpo_neoplasm_data,
        "commercial_panels": fetch_commercial_panels_data,
        "cosmic": fetch_cosmic_data,
        "somatic_commercial": fetch_somatic_commercial_data,
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
    if bed_config.get("germline", False) and "include_germline" in df.columns:
        bed_path = output_dir / "germline_panel.bed"
        create_bed_file(df, bed_path, "include_germline")

    if bed_config.get("somatic", False) and "include_somatic" in df.columns:
        bed_path = output_dir / "somatic_panel.bed"
        create_bed_file(df, bed_path, "include_somatic")

    if bed_config.get("combined", False) and "include_any" in df.columns:
        bed_path = output_dir / "combined_panel.bed"
        create_bed_file(df, bed_path, "include_any")


def display_summary(df: pd.DataFrame, config: dict[str, Any]) -> None:
    """Display pipeline summary statistics."""
    console.print("\n[bold blue]Pipeline Summary[/bold blue]")

    # Basic statistics
    table = Table(title="Gene Panel Statistics")
    table.add_column("Metric", style="cyan")
    table.add_column("Count", style="green")

    total_genes = len(df)
    table.add_row("Total genes", str(total_genes))

    if "include_germline" in df.columns:
        germline_count = df["include_germline"].sum()
        table.add_row("Germline panel genes", str(germline_count))

    if "include_somatic" in df.columns:
        somatic_count = df["include_somatic"].sum()
        table.add_row("Somatic panel genes", str(somatic_count))

    if "include_any" in df.columns:
        any_count = df["include_any"].sum()
        table.add_row("Total included genes", str(any_count))

    # Annotation statistics
    if "chromosome" in df.columns:
        annotated_count = (~df["chromosome"].isna()).sum()
        table.add_row("Genes with coordinates", str(annotated_count))

    if "mane_transcript" in df.columns:
        mane_count = (~df["mane_transcript"].isna()).sum()
        table.add_row("Genes with MANE transcripts", str(mane_count))

    console.print(table)

    # Top scoring genes
    if "total_score" in df.columns and total_genes > 0:
        console.print("\n[bold blue]Top 10 Scoring Genes[/bold blue]")
        top_genes = df.nlargest(10, "total_score")

        top_table = Table()
        top_table.add_column("Gene", style="cyan")
        top_table.add_column("Total Score", style="green")
        top_table.add_column("Germline", style="yellow")
        top_table.add_column("Somatic", style="yellow")

        for _, row in top_genes.iterrows():
            top_table.add_row(
                row["approved_symbol"],
                f"{row['total_score']:.2f}",
                f"{row.get('germline_score', 0):.2f}",
                f"{row.get('somatic_score', 0):.2f}",
            )

        console.print(top_table)


if __name__ == "__main__":
    app()
