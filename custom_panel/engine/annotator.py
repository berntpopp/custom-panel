"""
Gene annotation engine.

This module provides functionality to annotate genes with genomic information
using HGNC and Ensembl APIs.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import pandas as pd

from ..core.cache_manager import CacheManager
from ..core.config_manager import ConfigManager
from ..core.dataframe_utils import (
    safe_column_count,
)
from ..core.ensembl_client import EnsemblClient
from ..core.hgnc_client import HGNCClient

logger = logging.getLogger(__name__)


class GeneAnnotator:
    """Gene annotation engine using HGNC and Ensembl APIs."""

    def __init__(self, config: dict[str, Any] | None = None):
        """
        Initialize the gene annotator.

        Args:
            config: Configuration dictionary
        """
        self.config = config or {}

        # Initialize cache manager if caching is enabled
        perf_config = self.config.get("performance", {})
        cache_enabled = perf_config.get("enable_caching", True)
        cache_dir = perf_config.get("cache_dir", ".cache")
        cache_ttl = perf_config.get("cache_ttl", 2592000)  # 30 days default

        self.cache_manager = None
        if cache_enabled:
            self.cache_manager = CacheManager(
                cache_dir=cache_dir, cache_ttl=cache_ttl, enabled=True
            )
            logger.info(
                f"Cache enabled at {cache_dir} with TTL {cache_ttl/86400:.1f} days"
            )

        # Initialize API clients
        api_config = self.config.get("apis", {})

        hgnc_config = api_config.get("hgnc", {})
        self.hgnc_client = HGNCClient(
            timeout=hgnc_config.get("timeout", 30),
            max_retries=hgnc_config.get("max_retries", 3),
            retry_delay=hgnc_config.get("retry_delay", 1.0),
            cache_manager=self.cache_manager,
        )

        ensembl_config = api_config.get("ensembl", {})
        self.ensembl_client = EnsemblClient(
            timeout=ensembl_config.get("timeout", 60),
            max_retries=ensembl_config.get("max_retries", 3),
            retry_delay=ensembl_config.get("retry_delay", 1.0),
            transcript_batch_size=ensembl_config.get("transcript_batch_size", 50),
            cache_manager=self.cache_manager,
        )

        # Performance settings
        perf_config = self.config.get("performance", {})
        self.max_workers = perf_config.get("max_workers", 4)
        self.batch_size = perf_config.get(
            "batch_size", 200
        )  # Large batch size for coordinate lookup

        # Annotation settings
        annotation_config = self.config.get("annotation", {})
        self.include_coordinates = annotation_config.get("genomic_coordinates", True)
        self.include_transcripts = annotation_config.get("transcript_info", True)
        self.include_mane = annotation_config.get("mane_transcripts", True)
        self.include_descriptions = annotation_config.get("gene_descriptions", True)

        # Coverage calculation settings
        self.transcript_padding = annotation_config.get("transcript_padding", 25)
        self.gene_padding = annotation_config.get("gene_padding", 5000)

        config_manager = ConfigManager(self.config)
        self.species = config_manager.get_nested("general", "species", default="human")

    def annotate_genes(self, gene_df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate genes with genomic information.

        Args:
            gene_df: DataFrame with unique genes to annotate

        Returns:
            DataFrame with annotation columns added
        """
        if gene_df.empty:
            logger.warning("No genes to annotate")
            return gene_df

        # Get unique genes to annotate
        unique_genes = gene_df["approved_symbol"].dropna().unique().tolist()
        logger.info(f"Annotating {len(unique_genes)} unique genes")

        # Skip standardization - symbols were already standardized before merging
        # Just create identity mapping
        standardized_symbols = {symbol: symbol for symbol in unique_genes}

        # Step 1: Get genomic annotations
        annotations = self._get_gene_annotations(list(standardized_symbols.values()))

        # Store transcript data for potential exon BED file generation
        self.transcript_data = {
            symbol: annot
            for symbol, annot in annotations.items()
            if annot and "all_transcripts" in annot
        }

        # Step 2: Add annotations to the DataFrame
        annotated_df = self._add_annotations_to_dataframe(
            gene_df, standardized_symbols, annotations
        )

        logger.info(f"Successfully annotated {len(annotated_df)} gene records")
        return annotated_df

    def standardize_gene_symbols(
        self, gene_symbols: list[str]
    ) -> dict[str, dict[str, str | None]]:
        """
        Standardize gene symbols using parallel HGNC batch API calls.

        This method is public to allow standardization before merging.

        Args:
            gene_symbols: List of gene symbols to standardize

        Returns:
            Dictionary mapping original symbols to dict containing approved_symbol and hgnc_id
        """
        logger.info(
            f"Standardizing {len(gene_symbols)} gene symbols with parallel HGNC batch API"
        )

        # Split genes into batches
        batches = [
            gene_symbols[i : i + self.batch_size]
            for i in range(0, len(gene_symbols), self.batch_size)
        ]

        if len(batches) == 1:
            # Single batch - process directly
            return self._standardize_single_batch(batches[0])

        # Multiple batches - use parallel processing
        logger.info(
            f"Processing {len(batches)} HGNC batches in parallel (max_workers={self.max_workers})"
        )

        standardized = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all batch jobs
            future_to_batch = {
                executor.submit(self._standardize_single_batch, batch): batch
                for batch in batches
            }

            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch = future_to_batch[future]
                try:
                    batch_result = future.result()
                    standardized.update(batch_result)
                    logger.info(f"✓ Completed HGNC batch of {len(batch)} symbols")
                except Exception as e:
                    logger.error(
                        f"✗ HGNC batch processing failed for {len(batch)} symbols: {e}"
                    )
                    # Fallback to individual processing
                    batch_result = self._standardize_individual_parallel(batch)
                    standardized.update(batch_result)

        # Log standardization results
        changed_symbols = {
            k: v for k, v in standardized.items() if k != v["approved_symbol"]
        }
        if changed_symbols:
            logger.info(f"✨ Standardized {len(changed_symbols)} gene symbols")
            for original, info in changed_symbols.items():
                logger.debug(
                    f"  {original} -> {info['approved_symbol']} (HGNC ID: {info['hgnc_id']})"
                )

        return standardized

    def _standardize_single_batch(
        self, batch: list[str]
    ) -> dict[str, dict[str, str | None]]:
        """
        Standardize a single batch of gene symbols.

        Args:
            batch: List of gene symbols to standardize

        Returns:
            Dictionary mapping symbols to standardization results
        """
        try:
            # Use batch standardization for improved performance
            return self.hgnc_client.standardize_symbols(batch)
        except Exception as e:
            logger.warning(
                f"Batch standardization failed: {e}, falling back to individual requests"
            )
            return self._standardize_individual_parallel(batch)

    def _standardize_individual_parallel(
        self, symbols: list[str]
    ) -> dict[str, dict[str, str | None]]:
        """
        Standardize symbols individually in parallel when batch fails.

        Args:
            symbols: List of gene symbols to standardize

        Returns:
            Dictionary mapping symbols to standardization results
        """
        logger.info(f"Standardizing {len(symbols)} symbols individually in parallel")

        standardized = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit individual standardization jobs
            future_to_symbol = {
                executor.submit(self._standardize_single_symbol, symbol): symbol
                for symbol in symbols
            }

            # Collect results
            for future in as_completed(future_to_symbol):
                symbol = future_to_symbol[future]
                try:
                    result = future.result()
                    standardized[symbol] = result
                except Exception as e:
                    logger.warning(f"Failed to standardize {symbol}: {e}")
                    standardized[symbol] = {
                        "approved_symbol": symbol,
                        "hgnc_id": None,
                    }

        return standardized

    def _standardize_single_symbol(self, symbol: str) -> dict[str, str | None]:
        """
        Standardize a single gene symbol.

        Args:
            symbol: Gene symbol to standardize

        Returns:
            Dictionary with approved_symbol and hgnc_id
        """
        try:
            standardized_symbol = self.hgnc_client.standardize_symbol(symbol)
            # For individual lookups, get the gene info to have HGNC ID
            gene_info = self.hgnc_client.get_gene_info(standardized_symbol)
            hgnc_id = gene_info.get("hgnc_id") if gene_info else None
            return {
                "approved_symbol": standardized_symbol,
                "hgnc_id": hgnc_id,
            }
        except Exception:
            return {
                "approved_symbol": symbol,
                "hgnc_id": None,
            }

    def _get_gene_annotations(
        self, gene_symbols: list[str]
    ) -> dict[str, dict[str, Any]]:
        """
        Get genomic annotations for genes using optimized parallel batch API calls.

        Args:
            gene_symbols: List of standardized gene symbols

        Returns:
            Dictionary mapping gene symbols to annotation data
        """
        logger.info(
            f"Fetching genomic annotations for {len(gene_symbols)} genes using parallel batch API"
        )

        # Split genes into batches
        batches = [
            gene_symbols[i : i + self.batch_size]
            for i in range(0, len(gene_symbols), self.batch_size)
        ]

        if len(batches) == 1:
            # Single batch - no need for parallelization
            return self._process_single_batch(batches[0])

        # Multiple batches - use parallel processing
        logger.info(
            f"Processing {len(batches)} batches in parallel (max_workers={self.max_workers})"
        )

        annotations = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all batch jobs
            future_to_batch = {
                executor.submit(self._process_single_batch, batch): batch
                for batch in batches
            }

            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch = future_to_batch[future]
                try:
                    batch_annotations = future.result()
                    annotations.update(batch_annotations)
                    logger.info(f"✓ Completed batch of {len(batch)} genes")
                except Exception as e:
                    logger.error(
                        f"✗ Batch processing failed for {len(batch)} genes: {e}"
                    )
                    # Fallback to sequential processing for this batch
                    batch_annotations = self._process_single_batch_fallback(batch)
                    annotations.update(batch_annotations)

        return annotations

    def _process_single_batch(self, batch: list[str]) -> dict[str, dict[str, Any]]:
        """
        Process a single batch of genes with full error handling.

        Args:
            batch: List of gene symbols to process

        Returns:
            Dictionary mapping gene symbols to annotation data
        """
        annotations = {}

        try:
            # Try batch call with transcript expansion first (most complete)
            batch_data = self.ensembl_client.get_symbols_data_batch(
                batch, self.species, expand=True
            )

            for symbol in batch:
                gene_data = batch_data.get(symbol)
                if gene_data:
                    annotations[
                        symbol
                    ] = self._build_gene_annotation_from_expanded_data(
                        symbol, gene_data
                    )
                else:
                    annotations[symbol] = self._build_empty_annotation(symbol)

        except Exception as e:
            logger.warning(
                f"Batch API request failed: {e}, trying without transcript expansion"
            )
            # Try again without transcript expansion (faster, less complete)
            try:
                batch_coords = self.ensembl_client.get_symbols_data_batch(
                    batch, self.species, expand=False
                )
                for symbol in batch:
                    coords = batch_coords.get(symbol)
                    if coords:
                        annotations[symbol] = self._build_gene_annotation(
                            symbol, coords
                        )
                    else:
                        annotations[symbol] = self._build_empty_annotation(symbol)
            except Exception as e2:
                logger.warning(
                    f"Batch coordinate request also failed: {e2}, falling back to individual requests"
                )
                # Final fallback - parallel individual requests
                annotations = self._process_individual_genes_parallel(batch)

        return annotations

    def _process_single_batch_fallback(
        self, batch: list[str]
    ) -> dict[str, dict[str, Any]]:
        """
        Fallback processing for a single batch that failed in parallel execution.

        Args:
            batch: List of gene symbols to process

        Returns:
            Dictionary mapping gene symbols to annotation data
        """
        logger.info(f"Processing fallback for batch of {len(batch)} genes")
        annotations = {}

        # Use individual processing as final fallback
        for symbol in batch:
            annotations[symbol] = self._get_individual_annotation(symbol)

        return annotations

    def _process_individual_genes_parallel(
        self, gene_symbols: list[str]
    ) -> dict[str, dict[str, Any]]:
        """
        Process individual genes in parallel when batch requests fail.

        Args:
            gene_symbols: List of gene symbols to process individually

        Returns:
            Dictionary mapping gene symbols to annotation data
        """
        logger.info(f"Processing {len(gene_symbols)} genes individually in parallel")

        annotations = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit individual gene jobs
            future_to_symbol = {
                executor.submit(self._get_individual_annotation, symbol): symbol
                for symbol in gene_symbols
            }

            # Collect results
            for future in as_completed(future_to_symbol):
                symbol = future_to_symbol[future]
                try:
                    annotation = future.result()
                    annotations[symbol] = annotation
                except Exception as e:
                    logger.warning(f"Failed to get annotation for {symbol}: {e}")
                    annotations[symbol] = self._build_empty_annotation(symbol)

        return annotations

    def _get_individual_annotation(self, gene_symbol: str) -> dict[str, Any]:
        """Get annotation for a single gene."""
        try:
            coords = self.ensembl_client.get_gene_coordinates(gene_symbol, self.species)
            if coords:
                return self._build_gene_annotation(gene_symbol, coords)
            else:
                return self._build_empty_annotation(gene_symbol)
        except Exception as e:
            logger.warning(f"Failed to get annotation for {gene_symbol}: {e}")
            return self._build_empty_annotation(gene_symbol)

    def _build_gene_annotation_from_expanded_data(
        self, gene_symbol: str, gene_data: dict[str, Any]
    ) -> dict[str, Any]:
        """Build complete gene annotation from expanded batch data."""
        annotation = {
            "gene_id": gene_data.get("gene_id"),
            "chromosome": gene_data.get("chromosome"),
            "gene_start": gene_data.get("start"),
            "gene_end": gene_data.get("end"),
            "gene_strand": gene_data.get("strand"),
            "biotype": gene_data.get("biotype"),
            "gene_description": gene_data.get("description"),
            "gene_size": None,
            "gene_coverage_with_padding": None,
            "canonical_transcript": None,
            "mane_select_transcript": None,
            "mane_select_refseq": None,
            "mane_clinical_transcript": None,
            "mane_clinical_refseq": None,
            "canonical_transcript_coverage": None,
            "mane_select_coverage": None,
            "mane_clinical_coverage": None,
        }

        # Calculate gene size and coverage
        if gene_data.get("start") and gene_data.get("end"):
            annotation["gene_size"] = gene_data["end"] - gene_data["start"] + 1
            annotation[
                "gene_coverage_with_padding"
            ] = self.ensembl_client.calculate_gene_coverage(
                gene_data["start"], gene_data["end"], self.gene_padding
            )

        # Extract transcript information from expanded data
        if self.include_transcripts:
            canonical_info = gene_data.get("canonical_transcript")
            if canonical_info:
                annotation["canonical_transcript"] = canonical_info.get("transcript_id")

            # Calculate canonical transcript coverage if we have the full transcript data
            canonical_full = gene_data.get("canonical_transcript_full")
            if canonical_full:
                annotation[
                    "canonical_transcript_coverage"
                ] = self.ensembl_client.calculate_transcript_coverage(
                    canonical_full, self.transcript_padding
                )

        if self.include_mane:
            # MANE Select transcript
            mane_select = gene_data.get("mane_select")
            if mane_select:
                annotation["mane_select_transcript"] = mane_select.get("transcript_id")
                annotation["mane_select_refseq"] = mane_select.get("refseq_match")

            # Calculate MANE Select coverage
            mane_select_full = gene_data.get("mane_select_full")
            if mane_select_full:
                annotation[
                    "mane_select_coverage"
                ] = self.ensembl_client.calculate_transcript_coverage(
                    mane_select_full, self.transcript_padding
                )

            # MANE Plus Clinical transcript
            mane_clinical = gene_data.get("mane_clinical")
            if mane_clinical:
                annotation["mane_clinical_transcript"] = mane_clinical.get(
                    "transcript_id"
                )
                annotation["mane_clinical_refseq"] = mane_clinical.get("refseq_match")

            # Calculate MANE Clinical coverage
            mane_clinical_full = gene_data.get("mane_clinical_full")
            if mane_clinical_full:
                annotation[
                    "mane_clinical_coverage"
                ] = self.ensembl_client.calculate_transcript_coverage(
                    mane_clinical_full, self.transcript_padding
                )

        return annotation

    def _build_gene_annotation(
        self, gene_symbol: str, coords: dict[str, Any]
    ) -> dict[str, Any]:
        """Build complete gene annotation from coordinate data (fallback method)."""
        annotation = {
            "gene_id": coords.get("gene_id"),
            "chromosome": coords.get("chromosome"),
            "gene_start": coords.get("start"),
            "gene_end": coords.get("end"),
            "gene_strand": coords.get("strand"),
            "biotype": coords.get("biotype"),
            "gene_description": coords.get("description"),
            "gene_size": None,
            "gene_coverage_with_padding": None,
            "canonical_transcript": None,
            "mane_select_transcript": None,
            "mane_select_refseq": None,
            "mane_clinical_transcript": None,
            "mane_clinical_refseq": None,
            "canonical_transcript_coverage": None,
            "mane_select_coverage": None,
            "mane_clinical_coverage": None,
        }

        # Calculate gene size and coverage
        if coords.get("start") and coords.get("end"):
            annotation["gene_size"] = coords["end"] - coords["start"] + 1
            annotation[
                "gene_coverage_with_padding"
            ] = self.ensembl_client.calculate_gene_coverage(
                coords["start"], coords["end"], self.gene_padding
            )

        # For fallback, we can't get transcript info without additional API calls
        # This method is only used when batch API fails
        logger.debug(
            f"Using fallback annotation for {gene_symbol} - transcript info not available"
        )

        return annotation

    def _build_empty_annotation(self, gene_symbol: str) -> dict[str, Any]:
        """Build empty annotation for genes not found."""
        return {
            "gene_id": None,
            "chromosome": None,
            "gene_start": None,
            "gene_end": None,
            "gene_strand": None,
            "biotype": None,
            "gene_description": None,
            "gene_size": None,
            "gene_coverage_with_padding": None,
            "canonical_transcript": None,
            "mane_select_transcript": None,
            "mane_select_refseq": None,
            "mane_clinical_transcript": None,
            "mane_clinical_refseq": None,
            "canonical_transcript_coverage": None,
            "mane_select_coverage": None,
            "mane_clinical_coverage": None,
        }

    def _add_annotations_to_dataframe(
        self,
        df: pd.DataFrame,
        standardized_symbols: dict[str, str],
        annotations: dict[str, dict[str, Any]],
    ) -> pd.DataFrame:
        """
        Add annotations to the DataFrame.

        Args:
            df: Original DataFrame
            standardized_symbols: Mapping of original to standardized symbols
            annotations: Gene annotations

        Returns:
            DataFrame with annotation columns added
        """
        # Create a copy to avoid modifying the original
        annotated_df = df.copy()

        # Update approved symbols with standardized versions
        annotated_df["approved_symbol"] = annotated_df["approved_symbol"].map(
            lambda x: standardized_symbols.get(x, x)
        )

        # Add HGNC IDs if not already present
        if (
            "hgnc_id" not in annotated_df.columns
            or safe_column_count(annotated_df, "hgnc_id") == 0
        ):
            hgnc_ids = []
            for symbol in annotated_df["approved_symbol"]:
                hgnc_id = self.hgnc_client.symbol_to_hgnc_id(symbol)
                hgnc_ids.append(hgnc_id or "")
            annotated_df["hgnc_id"] = hgnc_ids

        # Add genomic annotations
        annotation_columns = [
            "gene_id",
            "chromosome",
            "gene_start",
            "gene_end",
            "gene_strand",
            "biotype",
            "gene_description",
            "gene_size",
            "gene_coverage_with_padding",
            "canonical_transcript",
            "mane_select_transcript",
            "mane_select_refseq",
            "mane_clinical_transcript",
            "mane_clinical_refseq",
            "canonical_transcript_coverage",
            "mane_select_coverage",
            "mane_clinical_coverage",
        ]

        for col in annotation_columns:
            values = []
            for symbol in annotated_df["approved_symbol"]:
                annotation = annotations.get(symbol, {})
                values.append(annotation.get(col))
            annotated_df[col] = values

        return annotated_df

    def get_annotation_summary(self, annotated_df: pd.DataFrame) -> dict[str, Any]:
        """
        Generate summary statistics for annotations.

        Args:
            annotated_df: Annotated DataFrame

        Returns:
            Summary statistics dictionary
        """
        if annotated_df.empty:
            return {"total_genes": 0}

        total_genes = len(annotated_df)

        summary: dict[str, Any] = {
            "total_genes": total_genes,
            "with_hgnc_id": safe_column_count(annotated_df, "hgnc_id"),
            "with_coordinates": safe_column_count(annotated_df, "chromosome"),
            "with_gene_id": safe_column_count(annotated_df, "gene_id"),
            "with_canonical_transcript": safe_column_count(
                annotated_df, "canonical_transcript"
            ),
            "with_mane_select": safe_column_count(
                annotated_df, "mane_select_transcript"
            ),
            "with_mane_clinical": safe_column_count(
                annotated_df, "mane_clinical_transcript"
            ),
            "with_description": safe_column_count(annotated_df, "gene_description"),
        }

        # Calculate percentages
        for key in [
            "with_hgnc_id",
            "with_coordinates",
            "with_gene_id",
            "with_canonical_transcript",
            "with_mane_select",
            "with_mane_clinical",
            "with_description",
        ]:
            percentage = (summary[key] / total_genes * 100) if total_genes > 0 else 0.0
            summary[f"{key}_percent"] = float(percentage)

        # Chromosome distribution
        if "chromosome" in annotated_df.columns:
            chrom_counts = annotated_df["chromosome"].value_counts().to_dict()
            summary["chromosome_distribution"] = dict(chrom_counts)

        # Biotype distribution
        if "biotype" in annotated_df.columns:
            biotype_counts = annotated_df["biotype"].value_counts().to_dict()
            summary["biotype_distribution"] = dict(biotype_counts)

        return summary

    def clear_caches(self) -> None:
        """Clear API client caches."""
        self.hgnc_client.clear_cache()
        self.ensembl_client.clear_cache()

        # Also clear disk cache
        if self.cache_manager:
            count = self.cache_manager.clear()
            logger.info(f"Cleared {count} disk cache entries")

        logger.info("Cleared annotation caches")

    def get_cache_stats(self) -> dict[str, Any]:
        """Get cache statistics from API clients."""
        stats = {
            "hgnc_cache": self.hgnc_client.get_cache_info(),
            "ensembl_cache": self.ensembl_client.get_cache_info(),
        }

        if self.cache_manager:
            stats["disk_cache"] = self.cache_manager.get_stats()

        return stats
