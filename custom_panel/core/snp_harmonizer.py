"""
SNP harmonization service for standardizing SNP identifiers and coordinates.

This module provides comprehensive SNP harmonization including:
- SNP ID standardization (prioritize rsIDs)
- Dual-build coordinate support (hg19 and hg38)
- Batch processing for efficiency
- Coordinate liftover using gnomAD (validation disabled for rsID processing)
"""

from __future__ import annotations

import logging
import time
from typing import Any, Optional

import pandas as pd
from tqdm import tqdm

from .ensembl_client import EnsemblClient
from .gnomad_client import GnomADClient

logger = logging.getLogger(__name__)


class HarmonizationMetrics:
    """Helper class for tracking and logging harmonization metrics."""

    def __init__(self) -> None:
        self.start_time = time.time()
        self.operation_times: dict[str, float] = {}
        self.error_details: list[dict[str, Any]] = []

    def log_operation_start(self, operation: str, snp_id: Optional[str] = None) -> None:
        """Log the start of an operation."""
        self.operation_times[operation] = time.time()
        if snp_id:
            logger.debug(f"▶️ Starting {operation} for SNP: {snp_id}")
        else:
            logger.debug(f"▶️ Starting {operation}")

    def log_operation_end(
        self,
        operation: str,
        success: bool,
        snp_id: Optional[str] = None,
        details: Optional[dict[str, Any]] = None,
    ) -> None:
        """Log the end of an operation with timing."""
        if operation in self.operation_times:
            duration = time.time() - self.operation_times[operation]
            status = "✅ SUCCESS" if success else "❌ FAILED"

            if snp_id:
                logger.debug(
                    f"{status} {operation} for SNP: {snp_id} in {duration:.3f}s",
                    extra={
                        "operation": operation,
                        "snp_id": snp_id,
                        "duration_seconds": duration,
                        "success": success,
                        "details": details or {},
                    },
                )
            else:
                logger.info(
                    f"{status} {operation} in {duration:.3f}s",
                    extra={
                        "operation": operation,
                        "duration_seconds": duration,
                        "success": success,
                        "details": details or {},
                    },
                )

            del self.operation_times[operation]

    def log_error(
        self,
        operation: str,
        error: Exception,
        snp_id: Optional[str] = None,
        context: Optional[dict[str, Any]] = None,
    ) -> None:
        """Log an error with context."""
        error_detail = {
            "operation": operation,
            "snp_id": snp_id,
            "error_type": type(error).__name__,
            "error_message": str(error),
            "context": context or {},
        }
        self.error_details.append(error_detail)

        logger.error(
            f"🚨 {operation} failed for SNP: {snp_id or 'Unknown'} - {error}",
            extra=error_detail,
        )

    def get_summary(self) -> dict[str, Any]:
        """Get a summary of harmonization metrics."""
        total_duration = time.time() - self.start_time
        return {
            "total_duration_seconds": total_duration,
            "error_count": len(self.error_details),
            "error_details": self.error_details,
        }


class SNPHarmonizationError(Exception):
    """Custom exception for SNP harmonization errors."""

    pass


class SNPHarmonizer:
    """
    Service for harmonizing SNP identifiers and coordinates across sources.

    Provides methods for:
    - Standardizing SNP IDs (prioritize rsIDs)
    - Ensuring dual-build coordinates
    - Batch processing for efficiency
    - Coordinate liftover using gnomAD (validation disabled for rsID processing)
    """

    def __init__(
        self,
        gnomad_client: GnomADClient,
        ensembl_client: EnsemblClient,
        config: Optional[dict[str, Any]] = None,
    ):
        """
        Initialize SNP harmonizer.

        Args:
            gnomad_client: Initialized gnomAD client
            ensembl_client: Initialized Ensembl client
            config: Optional configuration dictionary
        """
        self.gnomad_client = gnomad_client
        self.ensembl_client = ensembl_client
        self.config = config or {}

        # Configuration options
        self.prefer_rsid = self.config.get("prefer_rsid", True)
        self.fallback_to_coordinates = self.config.get("fallback_to_coordinates", True)
        self.batch_size = self.config.get("ensembl_batch_size", 25)

        # Statistics
        self.stats = {
            "total_processed": 0,
            "rsid_resolved": 0,
            "coordinates_lifted": 0,
            "errors": 0,
        }

        # Initialize metrics tracker
        self.metrics = HarmonizationMetrics()

        logger.info(
            "🔧 Initialized SNP harmonizer",
            extra={
                "prefer_rsid": self.prefer_rsid,
                "fallback_to_coordinates": self.fallback_to_coordinates,
                "validation_disabled": "gnomAD validation disabled for rsID processing",
            },
        )

    def harmonize_snp_batch(self, snp_df: pd.DataFrame) -> pd.DataFrame:
        """
        Harmonize a batch of SNPs using efficient batch processing.

        Args:
            snp_df: DataFrame with SNP data

        Returns:
            DataFrame with harmonized SNPs
        """
        if snp_df.empty:
            logger.warning("⚠️ Empty SNP DataFrame provided")
            return snp_df

        # Reset metrics for new batch
        self.metrics = HarmonizationMetrics()
        self.metrics.log_operation_start("SNP harmonization batch")

        logger.info(
            f"🔄 Harmonizing {len(snp_df)} SNPs",
            extra={
                "snp_count": len(snp_df),
                "config": {
                    "prefer_rsid": self.prefer_rsid,
                    "fallback_to_coordinates": self.fallback_to_coordinates,
                    "validation_disabled": "gnomAD validation disabled for rsID processing",
                },
            },
        )

        # Step 1: Collect all rsIDs that need coordinate resolution
        rsids_to_resolve = []
        for idx, row in snp_df.iterrows():
            snp_id = row.get("snp", f"index_{idx}")
            # Check if we need to resolve coordinates from rsID
            if (
                snp_id.startswith("rs")
                and snp_id != "rs"
                and not self._has_coordinates(row.to_dict())
            ):
                rsids_to_resolve.append(snp_id)

        # Step 2: Batch resolve coordinates for all rsIDs if any
        coordinate_cache = {}
        if rsids_to_resolve:
            logger.info(
                f"📡 Batch resolving coordinates for {len(rsids_to_resolve)} rsIDs using Ensembl API"
            )
            try:
                variations = self.ensembl_client.get_variations_batch(
                    rsids_to_resolve, batch_size=self.batch_size
                )

                for rsid in rsids_to_resolve:
                    if rsid in variations:
                        variation_data = variations[rsid]

                        # Extract coordinates for both builds
                        for build, assembly in [("hg38", "GRCh38"), ("hg19", "GRCh37")]:
                            coordinates = (
                                self.ensembl_client.extract_coordinates_from_variation(
                                    variation_data, preferred_assembly=assembly
                                )
                            )

                            if coordinates:
                                coordinate_cache[f"{rsid}_{build}"] = coordinates
                                logger.debug(
                                    f"✅ Cached coordinates for {rsid} ({assembly}): {coordinates.get('chromosome')}:{coordinates.get('start')}"
                                )

            except Exception as e:
                logger.error(f"❌ Failed to batch resolve coordinates: {e}")

        # Step 3: Process each SNP with pre-resolved coordinates
        harmonized_snps = []
        successful_harmonizations = 0

        for idx, row in tqdm(
            snp_df.iterrows(), total=len(snp_df), desc="Harmonizing SNPs"
        ):
            snp_id = row.get("snp", f"index_{idx}")
            try:
                harmonized_snp = self._harmonize_single_snp(row, coordinate_cache)
                if harmonized_snp is not None:
                    harmonized_snps.append(harmonized_snp)
                    if not harmonized_snp.get("harmonization_error"):
                        successful_harmonizations += 1
                self.stats["total_processed"] += 1

            except Exception as e:
                self.metrics.log_error(
                    "SNP harmonization",
                    e,
                    snp_id,
                    {"row_index": idx, "snp_data": row.to_dict()},
                )
                self.stats["errors"] += 1
                # Include original SNP with error flag
                error_snp = row.to_dict()
                error_snp["harmonization_error"] = str(e)
                harmonized_snps.append(error_snp)

        if harmonized_snps:
            result_df = pd.DataFrame(harmonized_snps)

            # Log comprehensive results
            self.metrics.log_operation_end(
                "SNP harmonization batch",
                True,
                details={
                    "total_snps": len(snp_df),
                    "successful_harmonizations": successful_harmonizations,
                    "errors": self.stats["errors"],
                    "success_rate": successful_harmonizations / len(snp_df)
                    if len(snp_df) > 0
                    else 0,
                },
            )

            logger.info(
                f"✅ Harmonized {len(result_df)} SNPs (success rate: {successful_harmonizations}/{len(snp_df)} = {successful_harmonizations/len(snp_df)*100:.1f}%)",
                extra={
                    "total_snps": len(snp_df),
                    "successful_harmonizations": successful_harmonizations,
                    "errors": self.stats["errors"],
                    "success_rate": successful_harmonizations / len(snp_df),
                    "harmonization_metrics": self.metrics.get_summary(),
                },
            )
            return result_df
        else:
            self.metrics.log_operation_end("SNP harmonization batch", False)
            logger.warning("⚠️ No SNPs were successfully harmonized")
            return pd.DataFrame()

    def _harmonize_single_snp(
        self, snp_row: pd.Series[Any], coordinate_cache: Optional[dict[str, Any]] = None
    ) -> Optional[dict[str, Any]]:
        """
        Harmonize a single SNP.

        Args:
            snp_row: Series with SNP data
            coordinate_cache: Pre-resolved coordinate cache

        Returns:
            Dictionary with harmonized SNP data
        """
        snp_info = snp_row.to_dict()

        # Step 1: Standardize SNP ID
        snp_info = self._standardize_snp_id(snp_info)

        # Step 2: Ensure dual coordinates
        snp_info = self._ensure_dual_coordinates(snp_info, coordinate_cache)

        # Step 3: Skip validation for rsID processing (only use gnomAD for coordinate liftover)
        # Validation is disabled for rsID processing as requested
        # gnomAD will only be used for coordinate liftover operations

        # Step 4: Add harmonization metadata
        snp_info["harmonized"] = True
        snp_info["harmonization_timestamp"] = pd.Timestamp.now().isoformat()

        return snp_info

    def _standardize_snp_id(self, snp_info: dict[str, Any]) -> dict[str, Any]:
        """
        Standardize SNP identifier, prioritizing rsIDs.

        Args:
            snp_info: SNP information dictionary

        Returns:
            Updated SNP information with standardized ID
        """
        current_id = snp_info.get("snp", "")

        # Check if we already have a proper rsID
        if current_id.startswith("rs") and current_id != "rs":
            # Already has rsID, keep it
            snp_info["rsid"] = current_id
            logger.debug(
                f"✅ SNP already has valid rsID: {current_id}",
                extra={"snp_id": current_id, "action": "keep_existing_rsid"},
            )
            return snp_info

        # Try to resolve rsID from coordinates
        if self.prefer_rsid and self._has_coordinates(snp_info):
            rsid = self._resolve_rsid_from_coordinates(snp_info)
            if rsid:
                snp_info["snp"] = rsid
                snp_info["rsid"] = rsid
                snp_info["original_id"] = current_id
                self.stats["rsid_resolved"] += 1
                logger.debug(
                    f"✅ Resolved rsID: {current_id} -> {rsid}",
                    extra={
                        "original_id": current_id,
                        "resolved_rsid": rsid,
                        "action": "rsid_resolved",
                    },
                )
                return snp_info

        # Fallback to coordinate-based ID if no rsID found
        if self.fallback_to_coordinates and self._has_coordinates(snp_info):
            coord_id = self._create_coordinate_id(snp_info)
            if coord_id:
                snp_info["snp"] = coord_id
                snp_info["coordinate_based_id"] = True
                snp_info["original_id"] = current_id
                logger.debug(
                    f"🔄 Using coordinate-based ID: {coord_id}",
                    extra={
                        "original_id": current_id,
                        "coordinate_id": coord_id,
                        "action": "coordinate_fallback",
                    },
                )
                return snp_info

        # Keep original ID if nothing else works
        logger.warning(
            f"⚠️ Could not standardize SNP ID: {current_id}",
            extra={"snp_id": current_id, "action": "failed_standardization"},
        )
        return snp_info

    def _ensure_dual_coordinates(
        self,
        snp_info: dict[str, Any],
        coordinate_cache: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        """
        Ensure SNP has both hg19 and hg38 coordinates using optimized resolution strategy.

        Optimization strategy:
        1. Use Ensembl batch variation API as primary method to get both builds
        2. Use gnomAD liftover only as fallback for missing specific builds
        3. Extract from coordinate-based SNP IDs as last resort

        Args:
            snp_info: SNP information dictionary
            coordinate_cache: Pre-resolved coordinate cache from batch processing

        Returns:
            Updated SNP information with dual coordinates
        """
        # Check what coordinates we have
        has_hg19 = self._has_build_coordinates(snp_info, "hg19")
        has_hg38 = self._has_build_coordinates(snp_info, "hg38")

        if has_hg19 and has_hg38:
            # Already have both, nothing to do
            return snp_info

        # Primary resolution: Try to resolve from rsID using Ensembl batch API first
        current_id = snp_info.get("snp", "")
        if current_id.startswith("rs") and current_id != "rs":
            coordinates_resolved = self._resolve_coordinates_from_rsid(
                snp_info, coordinate_cache
            )
            if coordinates_resolved:
                logger.debug(
                    f"✅ Resolved coordinates from rsID: {current_id}",
                    extra={"snp_id": current_id, "action": "coordinates_resolved"},
                )
                self.stats["coordinates_lifted"] += 1
                # Check if we now have both builds after resolution
                has_hg19 = self._has_build_coordinates(snp_info, "hg19")
                has_hg38 = self._has_build_coordinates(snp_info, "hg38")
                if has_hg19 and has_hg38:
                    return snp_info

        # Fallback: Use gnomAD liftover only for missing specific builds
        if has_hg38 and not has_hg19:
            # Only use liftover if Ensembl didn't provide hg19 coordinates
            self._liftover_coordinates(snp_info, "hg38", "hg19")
        elif has_hg19 and not has_hg38:
            # Only use liftover if Ensembl didn't provide hg38 coordinates
            self._liftover_coordinates(snp_info, "hg19", "hg38")

        # Last resort: Try to get coordinates from SNP ID if it's coordinate-based
        if not (
            self._has_build_coordinates(snp_info, "hg19")
            or self._has_build_coordinates(snp_info, "hg38")
        ):
            self._extract_coordinates_from_id(snp_info)

        return snp_info

    # NOTE: _validate_variant method is disabled for rsID processing
    # gnomAD validation is no longer used for rsID processing
    # gnomAD is only used for coordinate liftover operations
    #
    # def _validate_variant(self, snp_info: dict[str, Any]) -> dict[str, Any]:
    #     """
    #     Validate variant against gnomAD database.
    #
    #     This method is disabled for rsID processing as requested.
    #     gnomAD is only used for coordinate liftover operations.
    #     """
    #     pass

    def _resolve_rsid_from_coordinates(self, snp_info: dict[str, Any]) -> Optional[str]:
        """Resolve rsID from coordinates using gnomAD."""
        # Try hg38 first
        if self._has_build_coordinates(snp_info, "hg38"):
            try:
                position = snp_info.get("hg38_start", 0)
                if position == "" or position is None:
                    position = 0
                position = int(position)

                rsid = self.gnomad_client.resolve_rsid(
                    chromosome=str(snp_info.get("hg38_chromosome", "")),
                    position=position,
                    ref_allele=self._extract_ref_allele(snp_info, "hg38"),
                    alt_allele=self._extract_alt_allele(snp_info, "hg38"),
                    build="GRCh38",
                )
                if rsid:
                    return rsid
            except (ValueError, TypeError) as e:
                logger.debug(f"Error converting hg38_start to int: {e}")

        # Try hg19 as fallback
        if self._has_build_coordinates(snp_info, "hg19"):
            try:
                position = snp_info.get("hg19_start", 0)
                if position == "" or position is None:
                    position = 0
                position = int(position)

                rsid = self.gnomad_client.resolve_rsid(
                    chromosome=str(snp_info.get("hg19_chromosome", "")),
                    position=position,
                    ref_allele=self._extract_ref_allele(snp_info, "hg19"),
                    alt_allele=self._extract_alt_allele(snp_info, "hg19"),
                    build="GRCh37",
                )
                if rsid:
                    return rsid
            except (ValueError, TypeError) as e:
                logger.debug(f"Error converting hg19_start to int: {e}")

        return None

    def _resolve_coordinates_from_rsid(
        self,
        snp_info: dict[str, Any],
        coordinate_cache: Optional[dict[str, Any]] = None,
    ) -> Optional[bool]:
        """Resolve coordinates from rsID using cached coordinates or Ensembl batch variation API."""
        current_id = snp_info.get("snp", "")
        if not current_id.startswith("rs") or current_id == "rs":
            return None

        # First try to use the coordinate cache if available
        if coordinate_cache:
            for build, assembly in [("hg38", "GRCh38"), ("hg19", "GRCh37")]:
                cache_key = f"{current_id}_{build}"
                if cache_key in coordinate_cache:
                    coordinates = coordinate_cache[cache_key]

                    chromosome = coordinates.get("chromosome")
                    position = coordinates.get("start")
                    ref_allele = coordinates.get("ref_allele")
                    alt_allele = coordinates.get("alt_allele")

                    if all([chromosome, position, ref_allele, alt_allele]):
                        # Update SNP info with coordinates
                        snp_info[f"{build}_chromosome"] = str(chromosome)
                        snp_info[f"{build}_start"] = int(position) if position else 0
                        snp_info[f"{build}_end"] = int(position) if position else 0
                        snp_info[f"{build}_strand"] = 1
                        snp_info[
                            f"{build}_allele_string"
                        ] = f"{ref_allele}/{alt_allele}"

                        logger.debug(
                            f"✅ Used cached coordinates for {current_id} from {assembly}: {chromosome}:{position}:{ref_allele}:{alt_allele}",
                            extra={
                                "snp_id": current_id,
                                "build": assembly,
                                "chromosome": chromosome,
                                "position": position,
                                "ref_allele": ref_allele,
                                "alt_allele": alt_allele,
                                "source": "cache",
                            },
                        )
                        return True

        # Fallback to individual API call if not in cache
        try:
            # Use Ensembl batch variation API to get variant info
            variations = self.ensembl_client._get_variations_single_batch([current_id])

            if current_id in variations:
                variation_data = variations[current_id]

                # Try to extract coordinates for both builds
                for build, assembly in [("hg38", "GRCh38"), ("hg19", "GRCh37")]:
                    coordinates = (
                        self.ensembl_client.extract_coordinates_from_variation(
                            variation_data, preferred_assembly=assembly
                        )
                    )

                    if coordinates:
                        chromosome = coordinates.get("chromosome")
                        position = coordinates.get("start")
                        ref_allele = coordinates.get("ref_allele")
                        alt_allele = coordinates.get("alt_allele")

                        if all([chromosome, position, ref_allele, alt_allele]):
                            # Update SNP info with coordinates
                            snp_info[f"{build}_chromosome"] = str(chromosome)
                            snp_info[f"{build}_start"] = (
                                int(position) if position else 0
                            )
                            snp_info[f"{build}_end"] = int(position) if position else 0
                            snp_info[f"{build}_strand"] = 1
                            snp_info[
                                f"{build}_allele_string"
                            ] = f"{ref_allele}/{alt_allele}"

                            logger.debug(
                                f"✅ Resolved coordinates for {current_id} from {assembly}: {chromosome}:{position}:{ref_allele}:{alt_allele}",
                                extra={
                                    "snp_id": current_id,
                                    "build": assembly,
                                    "chromosome": chromosome,
                                    "position": position,
                                    "ref_allele": ref_allele,
                                    "alt_allele": alt_allele,
                                    "source": "ensembl_api",
                                },
                            )
                            return True

        except Exception as e:
            logger.debug(
                f"Failed to resolve coordinates from rsID {current_id} using Ensembl: {e}"
            )

        return None

    def _liftover_coordinates(
        self, snp_info: dict[str, Any], source_build: str, target_build: str
    ) -> bool:
        """Perform coordinate liftover using gnomAD."""
        # Get source coordinates
        source_chromosome = snp_info.get(f"{source_build}_chromosome")
        source_position = snp_info.get(f"{source_build}_start")
        source_ref = self._extract_ref_allele(snp_info, source_build)
        source_alt = self._extract_alt_allele(snp_info, source_build)

        # Validate and convert position
        if source_position == "" or source_position is None:
            logger.debug(f"Missing or empty position for {source_build}")
            return False

        try:
            source_position = int(source_position)
        except (ValueError, TypeError):
            logger.debug(
                f"Invalid position value for {source_build}: {source_position}"
            )
            return False

        if not all([source_chromosome, source_position, source_ref, source_alt]):
            logger.debug(f"Missing source coordinates for {source_build}")
            return False

        # Create source variant ID
        source_variant_id = (
            f"{source_chromosome}-{source_position}-{source_ref}-{source_alt}"
        )

        # Perform liftover
        build_map = {"hg19": "GRCh37", "hg38": "GRCh38"}
        source_build_name = build_map[source_build]
        target_build_name = build_map[target_build]

        result = self.gnomad_client.liftover_coordinates(
            source_variant_id=source_variant_id,
            source_build=source_build_name,
            target_build=target_build_name,
        )

        if result:
            (
                target_chromosome,
                target_variant_id,
                target_position,
                target_ref,
                target_alt,
            ) = result

            # Update SNP info with target coordinates
            snp_info[f"{target_build}_chromosome"] = target_chromosome
            snp_info[f"{target_build}_start"] = target_position
            snp_info[
                f"{target_build}_end"
            ] = target_position  # SNVs have same start/end
            snp_info[f"{target_build}_strand"] = 1  # Default positive strand
            snp_info[f"{target_build}_allele_string"] = f"{target_ref}/{target_alt}"

            self.stats["coordinates_lifted"] += 1
            logger.debug(
                f"Lifted coordinates: {source_variant_id} -> {target_variant_id}"
            )
            return True

        return False

    def _has_coordinates(self, snp_info: dict[str, Any]) -> bool:
        """Check if SNP has any coordinate information."""
        return (
            self._has_build_coordinates(snp_info, "hg19")
            or self._has_build_coordinates(snp_info, "hg38")
            or
            # Check for generic coordinate columns
            all(col in snp_info for col in ["chromosome", "position"])
        )

    def _has_build_coordinates(self, snp_info: dict[str, Any], build: str) -> bool:
        """Check if SNP has coordinates for specific build."""
        required_cols = [f"{build}_chromosome", f"{build}_start"]
        return all(
            col in snp_info and snp_info[col] is not None for col in required_cols
        )

    # NOTE: _validate_coordinates method is disabled for rsID processing
    # gnomAD validation is no longer used for rsID processing
    # gnomAD is only used for coordinate liftover operations
    #
    # def _validate_coordinates(self, snp_info: dict[str, Any], build: str) -> bool:
    #     """Validate coordinates against gnomAD."""
    #     # This method is disabled for rsID processing as requested
    #     # gnomAD is only used for coordinate liftover operations
    #     return True

    def _extract_ref_allele(self, snp_info: dict[str, Any], build: str) -> str:
        """Extract reference allele from SNP info."""
        # Try build-specific allele string first
        allele_string = snp_info.get(f"{build}_allele_string", "")
        if allele_string and "/" in allele_string:
            return allele_string.split("/")[0]

        # Try generic reference allele columns
        for col in ["ref_allele", "reference_allele", "ref"]:
            if col in snp_info and snp_info[col]:
                return str(snp_info[col])

        return ""

    def _extract_alt_allele(self, snp_info: dict[str, Any], build: str) -> str:
        """Extract alternate allele from SNP info."""
        # Try build-specific allele string first
        allele_string = snp_info.get(f"{build}_allele_string", "")
        if allele_string and "/" in allele_string:
            parts = allele_string.split("/")
            return parts[1] if len(parts) > 1 else ""

        # Try generic alternate allele columns
        for col in ["alt_allele", "alternate_allele", "alt", "effect_allele"]:
            if col in snp_info and snp_info[col]:
                return str(snp_info[col])

        return ""

    def _create_coordinate_id(self, snp_info: dict[str, Any]) -> Optional[str]:
        """Create coordinate-based ID from SNP info."""
        # Try hg38 first
        if self._has_build_coordinates(snp_info, "hg38"):
            chromosome = snp_info.get("hg38_chromosome")
            position = snp_info.get("hg38_start")
            ref = self._extract_ref_allele(snp_info, "hg38")
            alt = self._extract_alt_allele(snp_info, "hg38")

            if all([chromosome, position, ref, alt]):
                return f"{chromosome}:{position}:{ref}:{alt}"

        # Try hg19 as fallback
        if self._has_build_coordinates(snp_info, "hg19"):
            chromosome = snp_info.get("hg19_chromosome")
            position = snp_info.get("hg19_start")
            ref = self._extract_ref_allele(snp_info, "hg19")
            alt = self._extract_alt_allele(snp_info, "hg19")

            if all([chromosome, position, ref, alt]):
                return f"{chromosome}:{position}:{ref}:{alt}"

        return None

    def _extract_coordinates_from_id(self, snp_info: dict[str, Any]) -> None:
        """Extract coordinates from coordinate-based SNP ID."""
        snp_id = snp_info.get("snp", "")

        # Check if it's a coordinate-based ID (chr:pos:ref:alt)
        if ":" in snp_id:
            parts = snp_id.split(":")
            if len(parts) == 4:
                chromosome, position, ref, alt = parts
                try:
                    position = int(position)

                    # Assume it's hg38 by default (can be configured)
                    default_build = "hg38"
                    snp_info[f"{default_build}_chromosome"] = chromosome
                    snp_info[f"{default_build}_start"] = position
                    snp_info[f"{default_build}_end"] = position
                    snp_info[f"{default_build}_strand"] = 1
                    snp_info[f"{default_build}_allele_string"] = f"{ref}/{alt}"

                    logger.debug(f"Extracted coordinates from ID: {snp_id}")

                except ValueError:
                    logger.warning(f"Invalid position in coordinate ID: {snp_id}")

    def get_stats(self) -> dict[str, Any]:
        """Get harmonization statistics."""
        return self.stats.copy()

    def reset_stats(self) -> None:
        """Reset statistics."""
        self.stats = {
            "total_processed": 0,
            "rsid_resolved": 0,
            "coordinates_lifted": 0,
            "errors": 0,
        }
