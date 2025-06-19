"""
Improved command-line interface for the custom-panel tool.

This module provides the main CLI commands using Typer with better
modularization and use of the new utilities.
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
from rich.table import Table

from .core.config_manager import ConfigManager
from .core.output_generator import generate_outputs
from .engine.pipeline import Pipeline

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


def load_config(config_file: Optional[str] = None) -> dict[str, Any]:
    """Load configuration from file using improved config management."""
    try:
        # Load default config first
        default_config_path = Path(__file__).parent / "config" / "default_config.yml"
        config = _load_config_file(default_config_path, "default")

        # Apply overrides if specified
        if config_file:
            override_config = _load_config_file(Path(config_file), "override")
            config = _merge_configs(config, override_config)
        else:
            # Check for local overrides
            local_config_path = Path("config.local.yml")
            if local_config_path.exists():
                local_config = _load_config_file(local_config_path, "local")
                config = _merge_configs(config, local_config)

        _log_final_config_debug(config)
        return config

    except Exception as e:
        typer.echo(f"Error loading configuration: {e}", err=True)
        raise typer.Exit(1) from e


def _load_config_file(config_path: Path, config_type: str) -> dict[str, Any]:
    """Load a single configuration file."""
    if not config_path.exists():
        if config_type == "default":
            typer.echo(
                "Default configuration file not found. Installation may be corrupted.",
                err=True,
            )
            raise typer.Exit(1)
        return {}

    typer.echo(f"Loading {config_type} configuration from: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    _log_hpo_config_debug(config, config_type)
    return config or {}


def _log_hpo_config_debug(config: dict[str, Any], config_type: str) -> None:
    """Log HPO configuration for debugging."""
    config_manager = ConfigManager(config)
    hpo_config = config_manager.get_source_config("HPO_Neoplasm")
    typer.echo(
        f"DEBUG: {config_type.title()} HPO config - "
        f"URL: {hpo_config.get('omim_genemap2_url')}, "
        f"Path: {hpo_config.get('omim_genemap2_path')}"
    )


def _log_final_config_debug(config: dict[str, Any]) -> None:
    """Log final merged HPO configuration for debugging."""
    config_manager = ConfigManager(config)
    hpo_final = config_manager.get_source_config("HPO_Neoplasm")
    typer.echo(
        f"DEBUG: Final merged HPO config - "
        f"URL: {hpo_final.get('omim_genemap2_url')}, "
        f"Path: {hpo_final.get('omim_genemap2_path')}"
    )


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
    config_file: Optional[str] = typer.Option(
        None, "--config-file", "-c", help="Configuration file path"
    ),
    output_dir: Optional[str] = typer.Option(
        None, "--output-dir", "-o", help="Output directory"
    ),
    score_threshold: Optional[float] = typer.Option(
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

    # Load and configure
    config = load_config(config_file)
    config_manager = ConfigManager(config)

    # Override configuration with command line arguments
    config_manager.override_with_cli_args(
        output_dir=output_dir,
        score_threshold=score_threshold,
        save_intermediate=save_intermediate,
        intermediate_format=intermediate_format,
        log_to_file=log_to_file,
        structured_output=structured_output,
    )

    output_dir_path = Path(config_manager.get_output_dir())

    console.print("[bold green]Starting custom-panel pipeline[/bold green]")
    console.print(f"Output directory: {output_dir_path}")

    if dry_run:
        console.print(
            "[yellow]Running in dry-run mode - no files will be generated[/yellow]"
        )

    try:
        # Initialize and run the pipeline
        pipeline = Pipeline(config_manager.to_dict(), output_dir_path)
        annotated_df = pipeline.run(show_progress=True)

        # Generate outputs if not a dry run
        if not dry_run:
            _generate_pipeline_outputs(pipeline, annotated_df, config_manager)

        # Display summary
        display_summary(annotated_df, config_manager.to_dict())

        # Print success messages
        _print_completion_messages(pipeline, dry_run)

    except Exception as e:
        console.print(f"[red]Pipeline failed: {e}[/red]")
        if log_level.upper() == "DEBUG":
            import traceback

            console.print(traceback.format_exc())
        raise typer.Exit(1) from e


def _generate_pipeline_outputs(
    pipeline: Pipeline, annotated_df: pd.DataFrame, config_manager: ConfigManager
) -> None:
    """Generate all pipeline outputs."""
    final_output_dir = pipeline.output_manager.get_final_output_dir()
    generate_outputs(
        df=annotated_df,
        config=config_manager.to_dict(),
        output_dir=final_output_dir,
        transcript_data=pipeline.transcript_data,
    )

    # Save run summary
    run_summary = pipeline.get_run_summary()
    summary_path = pipeline.output_manager.run_dir / "run_summary.json"
    with open(summary_path, "w") as f:
        import json

        json.dump(run_summary, f, indent=2)

    # Cleanup old runs if using structured output
    pipeline.cleanup_old_runs()


def _print_completion_messages(pipeline: Pipeline, dry_run: bool) -> None:
    """Print appropriate completion messages."""
    if not dry_run:
        if pipeline.output_manager.use_structured:
            console.print(
                f"[bold green]Pipeline completed successfully! Results saved to: {pipeline.output_manager.run_dir}[/bold green]"
            )
            if pipeline.output_manager.intermediate_enabled:
                console.print("[blue]Intermediate files saved for debugging[/blue]")
            if pipeline.output_manager.file_logging_enabled:
                log_dir = (
                    pipeline.output_manager.run_dir
                    / pipeline.output_manager.subdirs.get("logs", "logs")
                )
                console.print(f"[blue]Logs saved to: {log_dir}[/blue]")
        else:
            console.print(
                f"[bold green]Pipeline completed successfully! Results saved to: {pipeline.output_manager.get_final_output_dir()}[/bold green]"
            )
    else:
        console.print("[bold yellow]Dry run completed successfully![/bold yellow]")


@app.command()
def fetch(
    source: str = typer.Argument(
        ...,
        help="Data source to fetch (use the pipeline instead - this command is deprecated)",
    ),
    config_file: Optional[str] = typer.Option(
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
    Fetch data from a single source (deprecated - use pipeline instead).
    """
    setup_logging(log_level)
    console.print(
        "[yellow]The fetch command is deprecated. Use 'run' command with specific source configuration instead.[/yellow]"
    )
    console.print(
        f"[yellow]Source '{source}' should be configured in your config file and run via the main pipeline.[/yellow]"
    )
    raise typer.Exit(1)


@app.command()
def config_check(
    config_file: Optional[str] = typer.Option(
        None, "--config-file", "-c", help="Configuration file path"
    ),
) -> None:
    """
    Validate and display configuration using improved config management.
    """
    config = load_config(config_file)
    config_manager = ConfigManager(config)

    console.print("[bold blue]Configuration Validation[/bold blue]")

    # Check data sources using config manager
    _display_data_sources_table(config_manager)

    # Check scoring configuration
    _display_scoring_configuration(config_manager)


def _display_data_sources_table(config_manager: ConfigManager) -> None:
    """Display data sources configuration table."""
    table = Table(title="Data Sources Configuration")
    table.add_column("Source", style="cyan")
    table.add_column("Enabled", style="green")
    table.add_column("Status", style="yellow")

    data_sources = config_manager.to_dict().get("data_sources", {})

    for source_name in data_sources.keys():
        enabled = config_manager.is_source_enabled(source_name)
        status = _get_source_status(config_manager, source_name, enabled)
        table.add_row(source_name, str(enabled), status)

    console.print(table)


def _get_source_status(
    config_manager: ConfigManager, source_name: str, enabled: bool
) -> str:
    """Get status string for a data source."""
    if not enabled:
        return "✗ Disabled"

    # Add specific validation for each source type
    if source_name == "inhouse_panels":
        return _validate_inhouse_panels_status(config_manager)
    elif source_name == "manual_curation":
        return _validate_manual_curation_status(config_manager)
    elif source_name == "hpo_neoplasm":
        return _validate_hpo_neoplasm_status(config_manager)
    else:
        return "✓ Configured"


def _validate_inhouse_panels_status(config_manager: ConfigManager) -> str:
    """Validate inhouse panels configuration."""
    source_config = config_manager.get_source_config("inhouse_panels")
    panels = source_config.get("panels", [])

    if not panels:
        return "⚠ No panels configured"

    missing_files = []
    for panel in panels:
        file_path = panel.get("file_path")
        if file_path and not Path(file_path).exists():
            missing_files.append(Path(file_path).name)

    if missing_files:
        return f"⚠ Missing files: {', '.join(missing_files)}"

    return "✓ Configured"


def _validate_manual_curation_status(config_manager: ConfigManager) -> str:
    """Validate manual curation configuration."""
    source_config = config_manager.get_source_config("manual_curation")
    lists = source_config.get("lists", [])

    if not lists:
        return "⚠ No lists configured"

    missing_files = []
    for list_item in lists:
        file_path = list_item.get("file_path")
        if file_path and not Path(file_path).exists():
            missing_files.append(Path(file_path).name)

    if missing_files:
        return f"⚠ Missing files: {', '.join(missing_files)}"

    return "✓ Configured"


def _validate_hpo_neoplasm_status(config_manager: ConfigManager) -> str:
    """Validate HPO neoplasm configuration."""
    source_config = config_manager.get_source_config("hpo_neoplasm")
    omim_file = source_config.get("omim_file_path")

    if omim_file and not Path(omim_file).exists():
        return f"⚠ OMIM file not found: {Path(omim_file).name}"

    specific_terms = source_config.get("specific_hpo_terms", [])
    if (
        not source_config.get("use_neoplasm_search", True)
        and not specific_terms
        and not omim_file
    ):
        return "⚠ No data sources configured"

    return "✓ Configured"


def _display_scoring_configuration(config_manager: ConfigManager) -> None:
    """Display scoring configuration."""
    scoring_config = config_manager.get_scoring_config()
    if scoring_config:
        console.print("\n[bold blue]Scoring Configuration[/bold blue]")

        score_threshold = config_manager.get_score_threshold()
        min_sources = config_manager.get_min_sources()

        console.print(f"Score threshold: {score_threshold or 'Not set'}")
        console.print(f"Minimum sources: {min_sources or 'Not set'}")


@app.command()
def search_panels(
    query: str = typer.Argument(..., help="Search term for panel names"),
    config_file: Optional[str] = typer.Option(
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


def display_summary(df: pd.DataFrame, config: dict[str, Any]) -> None:
    """Display pipeline summary statistics using improved utilities."""
    from .core.dataframe_utils import (
        safe_bool_count,
        safe_column_count,
    )

    console.print("\n[bold blue]Pipeline Summary[/bold blue]")

    # Basic statistics
    table = Table(title="Gene Panel Statistics")
    table.add_column("Metric", style="cyan")
    table.add_column("Count", style="green")

    total_genes = len(df)
    table.add_row("Total genes", str(total_genes))

    # Use improved utilities for column checks
    included_count = safe_bool_count(df, "include")
    table.add_row("Panel genes", str(included_count))

    annotated_count = safe_column_count(df, "chromosome")
    table.add_row("Genes with coordinates", str(annotated_count))

    mane_select_count = safe_column_count(df, "mane_select_transcript")
    mane_clinical_count = safe_column_count(df, "mane_clinical_transcript")
    table.add_row("Genes with MANE Select", str(mane_select_count))
    table.add_row("Genes with MANE Clinical", str(mane_clinical_count))

    console.print(table)

    # Top scoring genes
    _display_top_scoring_genes(df)


def _display_top_scoring_genes(df: pd.DataFrame) -> None:
    """Display top scoring genes table."""
    if "score" in df.columns and len(df) > 0:
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
