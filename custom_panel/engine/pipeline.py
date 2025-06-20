"""
Improved pipeline orchestrator for the custom-panel tool.

This module provides the main Pipeline class that encapsulates the entire
data processing workflow from fetching sources to annotation with better
modularization and use of DRY principles.
"""

import logging
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from ..core.config_manager import ConfigManager
from ..core.output_manager import OutputManager
from ..core.utils import normalize_count
from ..engine.annotator import GeneAnnotator
from ..engine.merger import PanelMerger
from ..engine.snp_annotator import SNPAnnotator
from ..sources.a_incidental_findings import fetch_acmg_incidental_data
from ..sources.b_manual_curation import fetch_manual_curation_data
from ..sources.g00_inhouse_panels import fetch_inhouse_panels_data
from ..sources.g01_panelapp import fetch_panelapp_data
from ..sources.g02_hpo import fetch_hpo_neoplasm_data
from ..sources.g03_commercial_panels import fetch_commercial_panels_data
from ..sources.g04_cosmic_germline import fetch_cosmic_germline_data
from ..sources.g05_clingen import fetch_clingen_data
from ..sources.g06_gencc import fetch_gencc_data
from ..sources_snp.ethnicity_snps import fetch_ethnicity_snps
from ..sources_snp.identity_snps import fetch_identity_snps
from ..sources_snp.manual_snps import fetch_manual_snps
from ..sources_snp.prs_snps import fetch_prs_snps

logger = logging.getLogger(__name__)
console = Console()


class Pipeline:
    """Orchestrates the gene panel curation pipeline with improved modularity."""

    def __init__(self, config: dict[str, Any], output_dir_path: Optional[Path] = None):
        """
        Initialize the pipeline with configuration.

        Args:
            config: Configuration dictionary
            output_dir_path: Optional output directory path override
        """
        self.config_manager = ConfigManager(config)
        self.output_dir_path = output_dir_path or Path(
            self.config_manager.get_output_dir()
        )
        self.output_manager = OutputManager(config, self.output_dir_path)
        self.annotator = GeneAnnotator(config)
        self.merger = PanelMerger(config)
        self.snp_annotator = SNPAnnotator(config)
        self.transcript_data: dict[str, Any] = {}
        self.snp_data: dict[str, pd.DataFrame] = {}

    def run(self, show_progress: bool = True) -> tuple[pd.DataFrame, dict[str, Any]]:
        """
        Execute the full gene panel curation pipeline.

        Args:
            show_progress: Whether to show progress indicators

        Returns:
            Tuple of (Final annotated DataFrame, transcript data)

        Raises:
            RuntimeError: If pipeline fails at any step
        """
        if show_progress:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console,
            ) as progress:
                return self._run_with_progress(progress)
        else:
            return self._run_without_progress()

    def _run_with_progress(
        self, progress: Progress
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Run pipeline with progress indicators."""
        task = progress.add_task("Starting pipeline...", total=None)

        # Step 1: Fetch data from all sources
        progress.update(task, description="Fetching data from sources...")
        raw_dataframes = self._fetch_all_sources()
        self._validate_fetched_data(raw_dataframes)
        progress.update(
            task, description=f"Fetched data from {len(raw_dataframes)} sources"
        )

        # Step 2: Centralized gene symbol standardization
        progress.update(task, description="Centralizing gene symbol standardization...")
        unified_df, symbol_map = self._standardize_symbols(
            raw_dataframes, progress, task
        )

        # Step 3: Pre-aggregate by source groups
        progress.update(task, description="Pre-aggregating source groups...")
        aggregated_sources = self._pre_aggregate_sources(unified_df, symbol_map)
        self._validate_aggregated_data(aggregated_sources)

        # Step 4: Merge and score pre-aggregated sources
        progress.update(
            task, description="Merging and scoring pre-aggregated sources..."
        )
        master_df = self.merger.create_master_list(
            aggregated_sources, self.output_manager
        )
        self._validate_master_data(master_df)

        # Step 5: Annotate the final master list
        progress.update(task, description="Annotating final gene list...")
        annotated_df = self._annotate_and_save(master_df)

        # Step 6: Process SNPs if enabled
        progress.update(task, description="Processing SNPs...")
        self._process_snps()

        progress.remove_task(task)
        return annotated_df, self.transcript_data

    def _run_without_progress(self) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Run pipeline without progress indicators."""
        # Step 1: Fetch data
        raw_dataframes = self._fetch_all_sources()
        self._validate_fetched_data(raw_dataframes)

        # Step 2: Standardize symbols
        unified_df, symbol_map = self._standardize_symbols(raw_dataframes)

        # Step 3: Pre-aggregate
        aggregated_sources = self._pre_aggregate_sources(unified_df, symbol_map)
        self._validate_aggregated_data(aggregated_sources)

        # Step 4: Merge and score
        master_df = self.merger.create_master_list(
            aggregated_sources, self.output_manager
        )
        self._validate_master_data(master_df)

        # Step 5: Annotate
        annotated_df = self._annotate_and_save(master_df)

        # Step 6: Process SNPs
        self._process_snps()

        return annotated_df, self.transcript_data

    def _fetch_all_sources(self) -> list[pd.DataFrame]:
        """Fetch data from all enabled sources."""
        # Define source functions
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

        dataframes = []
        for source_name, fetch_function in source_functions.items():
            if self.config_manager.is_source_enabled(source_name):
                dataframes.extend(
                    self._fetch_single_source(source_name, fetch_function)
                )
            else:
                console.print(
                    f"[yellow]Skipping disabled source: {source_name}[/yellow]"
                )

        return dataframes

    def _fetch_single_source(
        self, source_name: str, fetch_function: Any
    ) -> list[pd.DataFrame]:
        """Fetch data from a single source with error handling."""
        try:
            console.print(f"Fetching from {source_name}...")
            df = fetch_function(self.config_manager.to_dict())

            if not df.empty:
                # Save raw source data
                self.output_manager.save_source_data(df, source_name)
                console.print(f"[green]✓ {source_name}: {len(df)} records[/green]")
                return [df]
            else:
                console.print(f"[yellow]⚠ {source_name}: No data[/yellow]")
                return []

        except Exception as e:
            console.print(f"[red]✗ {source_name}: {e}[/red]")
            logger.exception(f"Error fetching from {source_name}")
            return []

    def _validate_fetched_data(self, dataframes: list[pd.DataFrame]) -> None:
        """Validate fetched data."""
        if not dataframes:
            raise RuntimeError("No data fetched from any source.")

    def _validate_aggregated_data(self, aggregated_sources: list[pd.DataFrame]) -> None:
        """Validate aggregated data."""
        if not aggregated_sources:
            raise RuntimeError("No valid genes after pre-aggregation.")

    def _validate_master_data(self, master_df: pd.DataFrame) -> None:
        """Validate master data."""
        if master_df.empty:
            raise RuntimeError("No genes in master list after merging.")

    def _standardize_symbols(
        self,
        raw_dataframes: list[pd.DataFrame],
        progress: Optional[Progress] = None,
        task_id: Any = None,
    ) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
        """Centralized gene symbol standardization."""
        # Combine all raw dataframes into one for efficient standardization
        unified_df = pd.concat(raw_dataframes, ignore_index=True)

        # Extract ALL unique symbols across all sources
        all_unique_symbols = unified_df["approved_symbol"].dropna().unique().tolist()
        logger.info(
            f"Total unique symbols across all sources: {len(all_unique_symbols)}"
        )

        # Make a single batch call to standardize ALL symbols
        if progress and task_id:
            progress.update(
                task_id,
                description=f"Standardizing {len(all_unique_symbols)} unique symbols...",
            )

        symbol_map = self.annotator.standardize_gene_symbols(all_unique_symbols)
        self._log_standardization_results(symbol_map, all_unique_symbols)

        # Apply standardization to unified dataframe
        unified_df = self._apply_symbol_standardization(unified_df, symbol_map)

        return unified_df, symbol_map

    def _log_standardization_results(
        self, symbol_map: dict[str, dict[str, Any]], all_symbols: list[str]
    ) -> None:
        """Log standardization results."""
        total_changed = sum(
            1
            for orig, info in symbol_map.items()
            if info["approved_symbol"] is not None and orig != info["approved_symbol"]
        )
        logger.info(f"Standardized {total_changed}/{len(all_symbols)} gene symbols")

    def _apply_symbol_standardization(
        self, unified_df: pd.DataFrame, symbol_map: dict[str, dict[str, Any]]
    ) -> pd.DataFrame:
        """Apply symbol standardization to the unified DataFrame."""
        # Create mapping dictionaries
        approved_symbol_map = {k: v["approved_symbol"] for k, v in symbol_map.items()}
        hgnc_id_map = {k: v["hgnc_id"] for k, v in symbol_map.items() if v["hgnc_id"]}

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

        return unified_df

    def _pre_aggregate_sources(
        self, unified_df: pd.DataFrame, symbol_map: dict[str, dict[str, Any]]
    ) -> list[pd.DataFrame]:
        """Pre-aggregate source groups."""
        # Add source_group column based on configuration
        source_to_group = self._build_source_group_mapping()
        unified_df["source_group"] = unified_df["source_name"].apply(
            lambda x: self._map_source_to_group(x, source_to_group)
        )

        # Pre-aggregate each source group
        aggregated_sources = []
        hgnc_id_map = {k: v["hgnc_id"] for k, v in symbol_map.items() if v["hgnc_id"]}

        for group_name, group_df in unified_df.groupby("source_group"):
            aggregated_df = self._aggregate_source_group(
                str(group_name), group_df, hgnc_id_map
            )
            aggregated_sources.append(aggregated_df)
            logger.info(
                f"Source group '{group_name}': {len(aggregated_df)} unique genes "
                f"from {group_df['source_name'].nunique()} sources"
            )

        return aggregated_sources

    def _build_source_group_mapping(self) -> dict[str, str]:
        """Build mapping from source names to groups."""
        source_to_group = {}

        for group_name in self.config_manager.to_dict().get("data_sources", {}):
            if self.config_manager.is_source_group(group_name):
                # This is a source group - map all its panels
                group_config = self.config_manager.get_source_config(group_name)
                for panel in group_config.get("panels", []):
                    if isinstance(panel, dict):
                        panel_name = panel.get("name")
                        if panel_name is not None:
                            source_to_group[panel_name] = group_name
            else:
                # Standalone source - it is its own group
                source_to_group[group_name] = group_name

        return source_to_group

    def _map_source_to_group(
        self, source_name: str, source_to_group: dict[str, str]
    ) -> str:
        """Map source name to its group."""
        # Try exact match first
        if source_name in source_to_group:
            return source_to_group[source_name]

        # Try prefix matching for sources with colons
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

    def _aggregate_source_group(
        self, group_name: str, group_df: pd.DataFrame, hgnc_id_map: dict[str, str]
    ) -> pd.DataFrame:
        """Aggregate genes within a source group."""
        gene_groups = []

        for gene_symbol, gene_df in group_df.groupby("approved_symbol"):
            aggregated_record = self._create_aggregated_record(
                str(gene_symbol), gene_df, group_name
            )
            gene_groups.append(aggregated_record)

        # Create DataFrame for this source group
        group_aggregated_df = pd.DataFrame(gene_groups)

        # Save intermediate aggregated data
        self.output_manager.save_standardized_data(
            group_aggregated_df,
            f"aggregated_{group_name}",
            {
                g: {"approved_symbol": g, "hgnc_id": hgnc_id_map.get(g)}
                for g in group_aggregated_df["approved_symbol"].unique()
            },
        )

        return group_aggregated_df

    def _create_aggregated_record(
        self, gene_symbol: str, gene_df: pd.DataFrame, group_name: str
    ) -> dict[str, Any]:
        """Create an aggregated record for a gene within a source group."""
        # Count internal sources
        internal_source_count = gene_df["source_name"].nunique()

        # Calculate internal confidence score if this is a source group
        if self.config_manager.is_source_group(group_name):
            normalization = self.config_manager.get_source_normalization(group_name)
            internal_confidence = normalize_count(
                internal_source_count,
                method=normalization.get("method", "linear"),
                params=normalization,
            )
        else:
            # Standalone sources have confidence of 1.0
            internal_confidence = 1.0

        # Create aggregated record for this gene in this source group
        return {
            "approved_symbol": gene_symbol,
            "hgnc_id": gene_df["hgnc_id"].iloc[0],  # Should be same for all
            "gene_name_reported": (
                gene_df["gene_name_reported"].iloc[0]
                if "gene_name_reported" in gene_df.columns
                else gene_symbol
            ),
            "source_name": group_name,  # Use group name as source for aggregated data
            "source_evidence_score": self.config_manager.get_source_evidence_score(
                group_name
            ),
            "source_details": f"{internal_source_count} sources in {group_name}",
            "source_group": group_name,
            "internal_source_count": internal_source_count,
            "internal_confidence_score": internal_confidence,
            "category": self.config_manager.get_source_category(group_name),
            "original_sources": [str(s) for s in gene_df["source_name"].unique()],
        }

    def _annotate_and_save(self, master_df: pd.DataFrame) -> pd.DataFrame:
        """Annotate the master DataFrame and save intermediate data."""
        annotated_df = self.annotator.annotate_genes(master_df)

        # Store transcript data for later use in output generation
        self.transcript_data = getattr(self.annotator, "transcript_data", {})

        # Save annotated data
        annotation_summary = self.annotator.get_annotation_summary(annotated_df)
        self.output_manager.save_annotated_data(annotated_df, annotation_summary)

        return annotated_df

    def _process_snps(self) -> None:
        """Process SNPs if enabled."""
        snp_config = self.config_manager.to_dict().get("snp_processing", {})

        if not snp_config.get("enabled", False):
            logger.info("SNP processing is disabled")
            return

        logger.info("Starting SNP processing...")

        # Define SNP source functions
        snp_source_functions = {
            "identity_and_ethnicity": fetch_identity_snps,
            "ethnicity": fetch_ethnicity_snps,
            "prs": fetch_prs_snps,
            "manual_snps": fetch_manual_snps,
        }

        # Fetch SNPs from all enabled sources
        for snp_type, fetch_function in snp_source_functions.items():
            snp_type_config = snp_config.get(snp_type, {})

            if snp_type_config.get("enabled", False):
                try:
                    console.print(f"Fetching {snp_type} SNPs...")
                    snp_df = fetch_function(self.config_manager.to_dict())

                    if snp_df is not None and not snp_df.empty:
                        # Annotate SNPs with genomic coordinates
                        annotated_snp_df = self.snp_annotator.annotate_snps(snp_df)

                        # Store SNP data
                        self.snp_data[snp_type] = annotated_snp_df

                        # Save SNP data using dedicated SNP save method
                        self.output_manager.save_snp_data(annotated_snp_df, snp_type)

                        console.print(
                            f"[green]✓ {snp_type}: {len(annotated_snp_df)} SNPs[/green]"
                        )
                        logger.info(
                            f"Successfully processed {len(annotated_snp_df)} {snp_type} SNPs"
                        )
                    else:
                        console.print(f"[yellow]⚠ {snp_type}: No SNPs found[/yellow]")

                except Exception as e:
                    console.print(f"[red]✗ {snp_type}: {e}[/red]")
                    logger.exception(f"Error processing {snp_type} SNPs")
            else:
                console.print(
                    f"[yellow]Skipping disabled SNP type: {snp_type}[/yellow]"
                )

        # Log SNP processing summary
        total_snps = sum(len(df) for df in self.snp_data.values())
        if total_snps > 0:
            logger.info(
                f"SNP processing completed: {total_snps} total SNPs across {len(self.snp_data)} categories"
            )
        else:
            logger.info("No SNPs were processed")

    def get_run_summary(self) -> dict[str, Any]:
        """Get summary of the pipeline run."""
        return self.output_manager.get_run_summary()

    def cleanup_old_runs(self) -> None:
        """Cleanup old pipeline runs if using structured output."""
        if self.config_manager.is_structured_output_enabled():
            self.output_manager.cleanup_old_runs()
