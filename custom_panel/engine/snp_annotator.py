"""
SNP annotation engine.

This module provides functionality to annotate SNPs with genomic information
using the Ensembl API. It handles batching and parallel processing for efficient
annotation of large SNP sets.
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import pandas as pd

from ..core.cache_manager import CacheManager
from ..core.config_manager import ConfigManager
from ..core.ensembl_client import EnsemblClient

logger = logging.getLogger(__name__)


class SNPAnnotator:
    """SNP annotation engine using Ensembl API for coordinate lookup."""

    def __init__(self, config: dict[str, Any] | None = None):
        """
        Initialize the SNP annotator.

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
                f"SNP Cache enabled at {cache_dir} with TTL {cache_ttl/86400:.1f} days"
            )

        # Initialize Ensembl client
        api_config = self.config.get("apis", {})
        ensembl_config = api_config.get("ensembl", {})
        self.ensembl_client = EnsemblClient(
            timeout=ensembl_config.get("timeout", 60),
            max_retries=ensembl_config.get("max_retries", 3),
            retry_delay=ensembl_config.get("retry_delay", 1.0),
            cache_manager=self.cache_manager,
        )

        # Performance settings
        perf_config = self.config.get("performance", {})
        self.max_workers = perf_config.get("max_workers", 4)
        self.batch_size = perf_config.get("batch_size", 100)  # Smaller batches for SNPs

        config_manager = ConfigManager(self.config)
        self.species = config_manager.get_nested("general", "species", default="human")

    def annotate_snps(self, snps_df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate SNPs DataFrame with genomic coordinates in parallel.
        
        This method replicates the R script pattern:
        ```r
        mutate(hg19 = snp_position_from_rs(snp, reference = "hg19"),
               hg38 = snp_position_from_rs(snp)) %>%
        unnest(c(hg19, hg38), names_sep = "_")
        ```

        Args:
            snps_df: DataFrame with SNPs (must have 'snp' column)

        Returns:
            DataFrame with SNP annotations including hg19 and hg38 coordinates
        """
        if snps_df.empty:
            logger.warning("No SNPs provided for annotation")
            return self._create_empty_annotation_dataframe()

        # Ensure we have the 'snp' column (matching R script)
        if "snp" not in snps_df.columns:
            if "rsid" in snps_df.columns:
                snps_df = snps_df.rename(columns={"rsid": "snp"})
            else:
                raise ValueError("DataFrame must contain 'snp' or 'rsid' column")

        # Get unique rsIDs to annotate
        unique_rsids = snps_df["snp"].dropna().unique().tolist()
        logger.info(f"Annotating {len(unique_rsids)} unique rsIDs with hg19 and hg38 coordinates")

        # Process rsIDs in batches using parallel execution for both assemblies
        hg38_annotations = self._get_snp_annotations_parallel(unique_rsids, assembly="GRCh38")
        hg19_annotations = self._get_snp_annotations_parallel(unique_rsids, assembly="GRCh37")

        # Add coordinate annotations to the original DataFrame
        annotated_df = self._add_coordinate_annotations_to_dataframe(
            snps_df, hg38_annotations, hg19_annotations
        )

        logger.info(
            f"Successfully annotated {len(annotated_df)} SNP records "
            f"({len(annotated_df[annotated_df['hg38_chromosome'].notna()])} with hg38, "
            f"{len(annotated_df[annotated_df['hg19_chromosome'].notna()])} with hg19)"
        )

        return annotated_df

    def _get_snp_annotations_parallel(
        self, rsids: list[str], assembly: str = "GRCh38"
    ) -> dict[str, dict[str, Any]]:
        """
        Get SNP annotations using parallel batch processing.

        Args:
            rsids: List of rsIDs to annotate
            assembly: Genome assembly (GRCh38 or GRCh37)

        Returns:
            Dictionary mapping rsIDs to annotation data
        """
        # Split rsIDs into batches
        batches = [
            rsids[i : i + self.batch_size]
            for i in range(0, len(rsids), self.batch_size)
        ]

        if len(batches) == 1:
            # Single batch - process directly
            return self._process_single_batch(batches[0])

        # Multiple batches - use parallel processing
        logger.info(
            f"Processing {len(batches)} SNP batches for {assembly} in parallel (max_workers={self.max_workers})"
        )

        annotations = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all batch jobs
            future_to_batch = {
                executor.submit(self._process_single_batch, batch, assembly): batch
                for batch in batches
            }

            # Collect results as they complete
            for future in as_completed(future_to_batch):
                batch = future_to_batch[future]
                try:
                    batch_annotations = future.result()
                    annotations.update(batch_annotations)
                    logger.info(f"✓ Completed SNP batch of {len(batch)} rsIDs")
                except Exception as e:
                    logger.error(
                        f"✗ SNP batch processing failed for {len(batch)} rsIDs: {e}"
                    )
                    # Fallback to individual processing
                    batch_annotations = self._process_individual_rsids(batch, assembly)
                    annotations.update(batch_annotations)

        return annotations

    def _process_single_batch(self, batch: list[str], assembly: str = "GRCh38") -> dict[str, dict[str, Any]]:
        """
        Process a single batch of rsIDs.

        Args:
            batch: List of rsIDs to process
            assembly: Genome assembly

        Returns:
            Dictionary mapping rsIDs to annotation data
        """
        annotations = {}

        # Process each rsID individually (Ensembl doesn't have batch SNP API)
        for rsid in batch:
            try:
                coords = self.ensembl_client.rsid_to_coordinates(rsid, self.species)
                if coords:
                    annotations[rsid] = self._build_snp_annotation(rsid, coords, assembly)
                else:
                    annotations[rsid] = self._build_empty_snp_annotation(rsid)
            except Exception as e:
                logger.warning(f"Failed to get coordinates for {rsid} ({assembly}): {e}")
                annotations[rsid] = self._build_empty_snp_annotation(rsid)

        return annotations

    def _process_individual_rsids(self, rsids: list[str], assembly: str = "GRCh38") -> dict[str, dict[str, Any]]:
        """
        Process individual rsIDs when batch processing fails.

        Args:
            rsids: List of rsIDs to process individually
            assembly: Genome assembly

        Returns:
            Dictionary mapping rsIDs to annotation data
        """
        logger.info(f"Processing {len(rsids)} rsIDs individually in parallel for {assembly}")

        annotations = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit individual rsID jobs
            future_to_rsid = {
                executor.submit(self._get_individual_snp_annotation, rsid, assembly): rsid
                for rsid in rsids
            }

            # Collect results
            for future in as_completed(future_to_rsid):
                rsid = future_to_rsid[future]
                try:
                    annotation = future.result()
                    annotations[rsid] = annotation
                except Exception as e:
                    logger.warning(f"Failed to get annotation for {rsid} ({assembly}): {e}")
                    annotations[rsid] = self._build_empty_snp_annotation(rsid)

        return annotations

    def _get_individual_snp_annotation(self, rsid: str, assembly: str = "GRCh38") -> dict[str, Any]:
        """Get annotation for a single rsID."""
        try:
            coords = self.ensembl_client.rsid_to_coordinates(rsid, self.species)
            if coords:
                return self._build_snp_annotation(rsid, coords, assembly)
            else:
                return self._build_empty_snp_annotation(rsid)
        except Exception as e:
            logger.warning(f"Failed to get annotation for {rsid} ({assembly}): {e}")
            return self._build_empty_snp_annotation(rsid)

    def _build_snp_annotation(
        self, rsid: str, coords: dict[str, Any], assembly: str = "GRCh38"
    ) -> dict[str, Any]:
        """
        Build SNP annotation from coordinate data.

        Args:
            rsid: RS identifier
            coords: Coordinate data from Ensembl
            assembly: Genome assembly

        Returns:
            Dictionary with SNP annotation data
        """
        return {
            "rsid": rsid,
            "chromosome": coords.get("chromosome"),
            "start": coords.get("start"),
            "end": coords.get("end"),
            "strand": coords.get("strand"),
            "allele_string": coords.get("allele_string"),
            "assembly": assembly,
        }

    def _build_empty_snp_annotation(self, rsid: str) -> dict[str, Any]:
        """
        Build empty annotation for SNPs not found.

        Args:
            rsid: RS identifier

        Returns:
            Dictionary with empty annotation data
        """
        return {
            "rsid": rsid,
            "chromosome": None,
            "start": None,
            "end": None,
            "strand": None,
            "allele_string": None,
            "assembly": None,
        }

    def _add_coordinate_annotations_to_dataframe(
        self, 
        snps_df: pd.DataFrame, 
        hg38_annotations: dict[str, dict[str, Any]], 
        hg19_annotations: dict[str, dict[str, Any]]
    ) -> pd.DataFrame:
        """
        Add coordinate annotations to SNPs DataFrame.
        
        This method replicates the R script pattern of adding hg19_ and hg38_ 
        prefixed columns for coordinates from both genome assemblies.

        Args:
            snps_df: Original SNPs DataFrame
            hg38_annotations: hg38 coordinate annotations
            hg19_annotations: hg19 coordinate annotations

        Returns:
            DataFrame with coordinate annotations added
        """
        result_df = snps_df.copy()
        
        # Add hg38 coordinates with hg38_ prefix
        for snp in result_df["snp"]:
            hg38_coords = hg38_annotations.get(snp, {})
            hg19_coords = hg19_annotations.get(snp, {})
            
            # Get the index for this SNP (handle potential duplicates)
            snp_indices = result_df[result_df["snp"] == snp].index
            
            for idx in snp_indices:
                # Add hg38 coordinates
                result_df.loc[idx, "hg38_chromosome"] = hg38_coords.get("chromosome")
                result_df.loc[idx, "hg38_start"] = hg38_coords.get("start")
                result_df.loc[idx, "hg38_end"] = hg38_coords.get("end")
                result_df.loc[idx, "hg38_strand"] = hg38_coords.get("strand")
                result_df.loc[idx, "hg38_allele_string"] = hg38_coords.get("allele_string")
                
                # Add hg19 coordinates
                result_df.loc[idx, "hg19_chromosome"] = hg19_coords.get("chromosome")
                result_df.loc[idx, "hg19_start"] = hg19_coords.get("start")
                result_df.loc[idx, "hg19_end"] = hg19_coords.get("end")
                result_df.loc[idx, "hg19_strand"] = hg19_coords.get("strand")
                result_df.loc[idx, "hg19_allele_string"] = hg19_coords.get("allele_string")
        
        return result_df

    def _create_annotation_dataframe(
        self, annotations: dict[str, dict[str, Any]]
    ) -> pd.DataFrame:
        """
        Create a DataFrame from SNP annotations.

        Args:
            annotations: Dictionary mapping rsIDs to annotation data

        Returns:
            DataFrame with SNP annotations
        """
        if not annotations:
            return self._create_empty_annotation_dataframe()

        annotation_records = list(annotations.values())
        df = pd.DataFrame(annotation_records)

        # Ensure consistent column order
        column_order = [
            "rsid",
            "chromosome", 
            "start",
            "end",
            "strand",
            "allele_string",
            "assembly",
        ]

        # Reorder columns if they exist
        existing_columns = [col for col in column_order if col in df.columns]
        other_columns = [col for col in df.columns if col not in column_order]
        df = df[existing_columns + other_columns]

        return df

    def _create_empty_annotation_dataframe(self) -> pd.DataFrame:
        """Create an empty annotation DataFrame with proper schema."""
        return pd.DataFrame(
            columns=[
                "snp",
                "source",
                "category",
                "hg38_chromosome",
                "hg38_start", 
                "hg38_end",
                "hg38_strand",
                "hg38_allele_string",
                "hg19_chromosome",
                "hg19_start", 
                "hg19_end",
                "hg19_strand",
                "hg19_allele_string",
            ]
        )

    def get_annotation_summary(self, annotated_df: pd.DataFrame) -> dict[str, Any]:
        """
        Generate summary statistics for SNP annotations.

        Args:
            annotated_df: Annotated DataFrame

        Returns:
            Summary statistics dictionary
        """
        if annotated_df.empty:
            return {"total_snps": 0}

        total_snps = len(annotated_df)
        with_coordinates = annotated_df["chromosome"].notna().sum()

        summary = {
            "total_snps": total_snps,
            "with_coordinates": int(with_coordinates),
            "with_coordinates_percent": (
                float(with_coordinates / total_snps * 100) if total_snps > 0 else 0.0
            ),
            "missing_coordinates": int(total_snps - with_coordinates),
        }

        # Chromosome distribution
        if "chromosome" in annotated_df.columns:
            chrom_counts = annotated_df["chromosome"].value_counts().to_dict()
            summary["chromosome_distribution"] = dict(chrom_counts)

        # Assembly information
        if "assembly" in annotated_df.columns:
            assembly_counts = annotated_df["assembly"].value_counts().to_dict()
            summary["assembly_distribution"] = dict(assembly_counts)

        return summary

    def clear_caches(self) -> None:
        """Clear API client caches."""
        self.ensembl_client.clear_cache()

        # Also clear disk cache
        if self.cache_manager:
            count = self.cache_manager.clear()
            logger.info(f"Cleared {count} SNP cache entries")

        logger.info("Cleared SNP annotation caches")

    def get_cache_stats(self) -> dict[str, Any]:
        """Get cache statistics."""
        stats = {
            "ensembl_cache": self.ensembl_client.get_cache_info(),
        }

        if self.cache_manager:
            stats["disk_cache"] = self.cache_manager.get_stats()

        return stats