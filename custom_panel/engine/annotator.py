"""
Gene annotation engine.

This module provides functionality to annotate genes with genomic information
using HGNC and Ensembl APIs.
"""

import logging
from typing import Any

import pandas as pd

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

        # Initialize API clients
        api_config = self.config.get("apis", {})

        hgnc_config = api_config.get("hgnc", {})
        self.hgnc_client = HGNCClient(
            timeout=hgnc_config.get("timeout", 30),
            max_retries=hgnc_config.get("max_retries", 3),
            retry_delay=hgnc_config.get("retry_delay", 1.0),
        )

        ensembl_config = api_config.get("ensembl", {})
        self.ensembl_client = EnsemblClient(
            timeout=ensembl_config.get("timeout", 30),
            max_retries=ensembl_config.get("max_retries", 3),
            retry_delay=ensembl_config.get("retry_delay", 1.0),
        )

        # Performance settings
        perf_config = self.config.get("performance", {})
        self.max_workers = perf_config.get("max_workers", 4)
        self.batch_size = perf_config.get("batch_size", 100)

        # Annotation settings
        annotation_config = self.config.get("annotation", {})
        self.include_coordinates = annotation_config.get("genomic_coordinates", True)
        self.include_transcripts = annotation_config.get("transcript_info", True)
        self.include_mane = annotation_config.get("mane_transcripts", True)
        self.include_descriptions = annotation_config.get("gene_descriptions", True)

        self.species = self.config.get("general", {}).get("species", "human")

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
        Standardize gene symbols using HGNC batch API.

        This method is public to allow standardization before merging.

        Args:
            gene_symbols: List of gene symbols to standardize

        Returns:
            Dictionary mapping original symbols to dict containing approved_symbol and hgnc_id
        """
        logger.info("Standardizing gene symbols with HGNC batch API")

        standardized = {}

        # Process genes in batches for better performance
        for i in range(0, len(gene_symbols), self.batch_size):
            batch = gene_symbols[i : i + self.batch_size]

            try:
                # Use batch standardization for improved performance
                batch_result = self.hgnc_client.standardize_symbols(batch)
                standardized.update(batch_result)

                # No need for additional individual lookups - the batch method already handles fallbacks internally
                # The batch method in standardize_symbols() already does individual lookups for symbols not found

            except Exception as e:
                logger.warning(
                    f"Batch standardization failed for batch: {e}, falling back to individual requests"
                )
                # Fall back to individual processing
                for symbol in batch:
                    try:
                        standardized_symbol = self.hgnc_client.standardize_symbol(
                            symbol
                        )
                        # For individual lookups, get the gene info to have HGNC ID
                        gene_info = self.hgnc_client.get_gene_info(standardized_symbol)
                        hgnc_id = gene_info.get("hgnc_id") if gene_info else None
                        standardized[symbol] = {
                            "approved_symbol": standardized_symbol,
                            "hgnc_id": hgnc_id,
                        }
                    except Exception as e2:
                        logger.warning(f"Failed to standardize {symbol}: {e2}")
                        standardized[symbol] = {
                            "approved_symbol": symbol,
                            "hgnc_id": None,
                        }

        # Log standardization results
        changed_symbols = {
            k: v for k, v in standardized.items() if k != v["approved_symbol"]
        }
        if changed_symbols:
            logger.info(f"Standardized {len(changed_symbols)} gene symbols")
            for original, info in changed_symbols.items():
                logger.debug(
                    f"  {original} -> {info['approved_symbol']} (HGNC ID: {info['hgnc_id']})"
                )

        return standardized

    def _get_gene_annotations(
        self, gene_symbols: list[str]
    ) -> dict[str, dict[str, Any]]:
        """
        Get genomic annotations for genes.

        Args:
            gene_symbols: List of standardized gene symbols

        Returns:
            Dictionary mapping gene symbols to annotation data
        """
        logger.info("Fetching genomic annotations from Ensembl")

        annotations = {}

        # Process genes in batches
        for i in range(0, len(gene_symbols), self.batch_size):
            batch = gene_symbols[i : i + self.batch_size]

            # Get batch coordinates (Ensembl supports batch requests)
            if len(batch) > 1:
                try:
                    batch_coords = self.ensembl_client.get_genes_coordinates(
                        batch, self.species
                    )
                    for symbol in batch:
                        coords = batch_coords.get(symbol)
                        if coords:
                            annotations[symbol] = self._build_gene_annotation(
                                symbol, coords
                            )
                        else:
                            annotations[symbol] = self._build_empty_annotation(symbol)
                except Exception as e:
                    logger.warning(
                        f"Batch coordinate request failed: {e}, falling back to individual requests"
                    )
                    # Fall back to individual requests
                    for symbol in batch:
                        annotations[symbol] = self._get_individual_annotation(symbol)
            else:
                # Single gene
                annotations[batch[0]] = self._get_individual_annotation(batch[0])

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

    def _build_gene_annotation(
        self, gene_symbol: str, coords: dict[str, Any]
    ) -> dict[str, Any]:
        """Build complete gene annotation from coordinate data."""
        annotation = {
            "gene_id": coords.get("gene_id"),
            "chromosome": coords.get("chromosome"),
            "gene_start": coords.get("start"),
            "gene_end": coords.get("end"),
            "gene_strand": coords.get("strand"),
            "biotype": coords.get("biotype"),
            "gene_description": coords.get("description"),
            "gene_size": None,
            "canonical_transcript": None,
            "mane_transcript": None,
            "mane_type": None,
        }

        # Calculate gene size
        if coords.get("start") and coords.get("end"):
            annotation["gene_size"] = coords["end"] - coords["start"] + 1

        # Get additional annotations if enabled
        gene_id = coords.get("gene_id")

        if self.include_transcripts and gene_id:
            try:
                canonical_info = self.ensembl_client.get_canonical_transcript(
                    gene_id, self.species
                )
                if canonical_info:
                    annotation["canonical_transcript"] = canonical_info.get(
                        "transcript_id"
                    )
            except Exception as e:
                logger.debug(
                    f"Failed to get canonical transcript for {gene_symbol}: {e}"
                )

        if self.include_mane:
            try:
                mane_info = self.ensembl_client.get_mane_transcript(
                    gene_symbol, self.species
                )
                if mane_info:
                    annotation["mane_transcript"] = mane_info.get("transcript_id")
                    annotation["mane_type"] = mane_info.get("mane_type")
            except Exception as e:
                logger.debug(f"Failed to get MANE transcript for {gene_symbol}: {e}")

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
            "canonical_transcript": None,
            "mane_transcript": None,
            "mane_type": None,
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
            or annotated_df["hgnc_id"].isna().all()
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
            "canonical_transcript",
            "mane_transcript",
            "mane_type",
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
            "with_hgnc_id": (
                ~annotated_df["hgnc_id"].isna() & (annotated_df["hgnc_id"] != "")
            ).sum(),
            "with_coordinates": (~annotated_df["chromosome"].isna()).sum(),
            "with_gene_id": (~annotated_df["gene_id"].isna()).sum(),
            "with_canonical_transcript": (
                ~annotated_df["canonical_transcript"].isna()
            ).sum(),
            "with_mane_transcript": (~annotated_df["mane_transcript"].isna()).sum(),
            "with_description": (~annotated_df["gene_description"].isna()).sum(),
        }

        # Calculate percentages
        for key in [
            "with_hgnc_id",
            "with_coordinates",
            "with_gene_id",
            "with_canonical_transcript",
            "with_mane_transcript",
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
        logger.info("Cleared annotation caches")

    def get_cache_stats(self) -> dict[str, Any]:
        """Get cache statistics from API clients."""
        return {
            "hgnc_cache": self.hgnc_client.get_cache_info(),
            "ensembl_cache": self.ensembl_client.get_cache_info(),
        }
