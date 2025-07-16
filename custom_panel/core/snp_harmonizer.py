"""
SNP harmonization service for standardizing SNP identifiers and coordinates.

This module provides simple, efficient SNP harmonization including:
- SNP ID standardization (prioritize rsIDs)
- hg38 coordinate resolution via Ensembl batch API
- Batch processing for efficiency
"""

from __future__ import annotations

import logging
from typing import Any, Optional

import pandas as pd

from .ensembl_client import EnsemblClient

logger = logging.getLogger(__name__)


class SNPHarmonizer:
    """
    Simple, elegant SNP harmonizer for hg38 coordinate resolution.

    Uses Ensembl batch API to resolve coordinates for rsIDs efficiently.
    """

    def __init__(
        self, ensembl_client: EnsemblClient, config: Optional[dict[str, Any]] = None
    ):
        """
        Initialize SNP harmonizer.

        Args:
            ensembl_client: Initialized Ensembl client
            config: Optional configuration dictionary
        """
        self.ensembl_client = ensembl_client
        self.config = config or {}
        self.batch_size = self.config.get("ensembl_batch_size", 25)

        # Statistics
        self.stats = {
            "total_processed": 0,
            "coordinates_resolved": 0,
            "errors": 0,
        }

        logger.info("🔧 Initialized simplified SNP harmonizer (hg38 only)")

    def harmonize_snp_batch(self, snp_df: pd.DataFrame) -> pd.DataFrame:
        """
        Harmonize a batch of SNPs using Ensembl batch API for hg38 coordinates.

        Args:
            snp_df: DataFrame with SNP data (must have 'snp' column with rsIDs)

        Returns:
            DataFrame with hg38 coordinates added
        """
        if snp_df.empty:
            logger.warning("⚠️ Empty SNP DataFrame provided")
            return snp_df

        logger.info(f"🔄 Harmonizing {len(snp_df)} SNPs to hg38 coordinates")

        # Step 1: Extract unique rsIDs that need coordinate resolution
        rsids_to_resolve = []
        for _, row in snp_df.iterrows():
            snp_id = row.get("snp", "")
            if snp_id.startswith("rs") and snp_id != "rs":
                rsids_to_resolve.append(snp_id)

        if not rsids_to_resolve:
            logger.warning("⚠️ No valid rsIDs found for coordinate resolution")
            return snp_df

        unique_rsids = list(set(rsids_to_resolve))
        logger.info(f"📡 Resolving coordinates for {len(unique_rsids)} unique rsIDs")

        # Step 2: Batch resolve coordinates using Ensembl API
        coordinate_cache = {}
        try:
            variations = self.ensembl_client.get_variations_batch(
                unique_rsids, batch_size=self.batch_size
            )

            for rsid in unique_rsids:
                if rsid in variations:
                    variation_data = variations[rsid]
                    coordinates = (
                        self.ensembl_client.extract_coordinates_from_variation(
                            variation_data, preferred_assembly="GRCh38"
                        )
                    )

                    if coordinates:
                        coordinate_cache[rsid] = coordinates
                        self.stats["coordinates_resolved"] += 1

        except Exception as e:
            logger.error(f"❌ Failed to batch resolve coordinates: {e}")

        # Step 3: Add coordinates to SNP data
        result_df = snp_df.copy()

        for idx, row in result_df.iterrows():
            snp_id = row.get("snp", "")
            if snp_id in coordinate_cache:
                coords = coordinate_cache[snp_id]
                # Use .at for scalar assignment to avoid mypy indexing issues
                result_df.at[idx, "hg38_chromosome"] = coords.get("chromosome")
                result_df.at[idx, "hg38_start"] = coords.get("start")
                result_df.at[idx, "hg38_end"] = coords.get("end")
                result_df.at[idx, "hg38_strand"] = coords.get("strand")
                result_df.at[idx, "hg38_allele_string"] = coords.get("allele_string")

        self.stats["total_processed"] = len(snp_df)

        logger.info(
            f"✅ Successfully resolved coordinates for {self.stats['coordinates_resolved']}/{len(unique_rsids)} rsIDs"
        )

        return result_df

    def get_stats(self) -> dict[str, Any]:
        """Get harmonization statistics."""
        return self.stats.copy()

    def reset_stats(self) -> None:
        """Reset statistics."""
        self.stats = {
            "total_processed": 0,
            "coordinates_resolved": 0,
            "errors": 0,
        }
