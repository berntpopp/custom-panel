"""
SNP harmonization service for standardizing SNP identifiers and coordinates.

This module provides simple, efficient SNP harmonization including:
- SNP ID standardization (prioritize rsIDs)
- hg38 coordinate resolution via Ensembl batch API
- Batch processing for efficiency
"""

from __future__ import annotations

import logging
import re
from typing import Any

import pandas as pd

from .ensembl_client import EnsemblClient

logger = logging.getLogger(__name__)


class SNPHarmonizer:
    """
    Simple, elegant SNP harmonizer for hg38 coordinate resolution.

    Uses Ensembl batch API to resolve coordinates for rsIDs efficiently.
    """

    def __init__(
        self,
        ensembl_client: EnsemblClient,
        config: dict[str, Any] | None = None,
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
                unique_rsids,
                batch_size=self.batch_size,
            )

            # Build mapping from original rsIDs to canonical rsIDs (handles merged rsIDs)
            rsid_mapping = {}
            for original_rsid in unique_rsids:
                # First check if original rsID is directly in response (not merged)
                if original_rsid in variations:
                    rsid_mapping[original_rsid] = original_rsid
                else:
                    # Check if it's a synonym of any canonical rsID
                    for canonical_rsid, variation_data in variations.items():
                        synonyms = variation_data.get("synonyms", [])
                        if original_rsid in synonyms:
                            rsid_mapping[original_rsid] = canonical_rsid
                            logger.debug(
                                f"📍 Mapped merged rsID {original_rsid} → {canonical_rsid}",
                            )
                            break

            for rsid in unique_rsids:
                canonical_rsid: str | None = rsid_mapping.get(rsid)
                if canonical_rsid and canonical_rsid in variations:
                    variation_data = variations[canonical_rsid]
                    coordinates = (
                        self.ensembl_client.extract_coordinates_from_variation(
                            variation_data,
                            preferred_assembly="GRCh38",
                        )
                    )

                    if coordinates:
                        coordinate_cache[rsid] = coordinates
                        self.stats["coordinates_resolved"] += 1

            # Fallback: Use source coordinates for unresolved rsIDs
            unresolved_rsids = [
                rsid for rsid in unique_rsids if rsid not in coordinate_cache
            ]
            if unresolved_rsids:
                source_coordinates = self._extract_source_coordinates(
                    unresolved_rsids,
                    snp_df,
                )
                for rsid, coords in source_coordinates.items():
                    if coords:
                        coordinate_cache[rsid] = coords
                        self.stats["coordinates_resolved"] += 1
                        logger.debug(f"📍 Used source coordinates for {rsid}")

        except Exception as e:
            logger.error(f"❌ Failed to batch resolve coordinates: {e}")

        # Step 3: Add coordinates to SNP data and generate VCF format IDs
        result_df = snp_df.copy()

        # Preserve original rsIDs in rsid column if not already present
        if "rsid" not in result_df.columns:
            result_df["rsid"] = result_df["snp"].apply(
                lambda x: x if str(x).startswith("rs") else pd.NA,
            )

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

                # Generate VCF format ID (chr:pos:ref:alt) when coordinates are available
                chromosome = coords.get("chromosome")
                start = coords.get("start")
                allele_string = coords.get("allele_string")

                if chromosome and start and allele_string:
                    # Parse allele string (format: "ref/alt" or "ref/alt1/alt2")
                    alleles = allele_string.split("/")
                    if len(alleles) >= 2:
                        ref_allele = alleles[0]
                        alt_allele = alleles[1]  # Use first alt allele
                        vcf_id = f"{chromosome}:{start}:{ref_allele}:{alt_allele}"
                        result_df.at[idx, "snp"] = vcf_id
                        logger.debug(f"Generated VCF ID: {vcf_id} for rsID: {snp_id}")

                # Ensure rsid column contains the original rsID
                if str(snp_id).startswith("rs"):
                    result_df.at[idx, "rsid"] = snp_id

        self.stats["total_processed"] = len(snp_df)

        logger.info(
            f"✅ Successfully resolved coordinates for {self.stats['coordinates_resolved']}/{len(unique_rsids)} rsIDs",
        )

        return result_df

    def _extract_source_coordinates(
        self,
        rsids: list[str],
        snp_df: pd.DataFrame,
    ) -> dict[str, dict[str, Any] | None]:
        """
        Extract genomic coordinates from source data when Ensembl lookup fails.

        Supports multiple coordinate formats:
        - NCBI RefSeq: NC_000010.11:94949144 (PharmGKB format)
        - Direct coordinates: chr10:94949144

        Args:
            rsids: List of rsIDs to find coordinates for
            snp_df: Original DataFrame with source coordinate information

        Returns:
            Dictionary mapping rsID to coordinate information
        """
        coordinate_cache: dict[str, dict[str, Any] | None] = {}

        for rsid in rsids:
            # Find matching row(s) in source data
            matching_rows = snp_df[snp_df["snp"] == rsid]
            if matching_rows.empty:
                coordinate_cache[rsid] = None
                continue

            row = matching_rows.iloc[0]  # Use first match

            # Look for coordinate information in various columns
            location = None
            for col in ["location", "genomic_location", "coordinates", "position"]:
                if col in row and pd.notna(row[col]) and str(row[col]).strip():
                    location = str(row[col]).strip()
                    break

            if not location:
                coordinate_cache[rsid] = None
                continue

            # Parse coordinate formats
            coords = self._parse_genomic_location(location)
            coordinate_cache[rsid] = coords

        return coordinate_cache

    def _parse_genomic_location(self, location: str) -> dict[str, Any] | None:
        """
        Parse genomic location from various formats.

        Supported formats:
        - NCBI RefSeq: NC_000010.11:94949144
        - Direct: chr10:94949144
        - With range: chr10:94949144-94949144

        Args:
            location: Genomic location string

        Returns:
            Dictionary with parsed coordinates or None if parsing fails
        """
        if not location:
            return None

        try:
            # Handle NCBI RefSeq format (NC_000010.11:94949144)
            if location.startswith("NC_"):
                # Extract chromosome number from NC_000010.11
                chr_match = re.match(r"NC_0+(\d+)\..*?:(\d+)", location)
                if chr_match:
                    chromosome = chr_match.group(1)
                    position = int(chr_match.group(2))

                    return {
                        "chromosome": chromosome,
                        "start": position,
                        "end": position,
                        "strand": 1,
                        "allele_string": "N/N",  # Unknown alleles
                        "assembly": "GRCh38",
                    }

            # Handle direct chromosome format (chr10:94949144)
            elif ":" in location:
                parts = location.split(":")
                if len(parts) >= 2:
                    chr_part = parts[0].replace("chr", "")
                    pos_part = parts[1]

                    # Handle position ranges (94949144-94949144)
                    if "-" in pos_part:
                        pos_start, pos_end = pos_part.split("-", 1)
                        start_pos = int(pos_start)
                        end_pos = int(pos_end)
                    else:
                        start_pos = end_pos = int(pos_part)

                    return {
                        "chromosome": chr_part,
                        "start": start_pos,
                        "end": end_pos,
                        "strand": 1,
                        "allele_string": "N/N",  # Unknown alleles
                        "assembly": "GRCh38",
                    }

        except (ValueError, AttributeError) as e:
            logger.debug(f"Failed to parse genomic location '{location}': {e}")

        return None

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
