"""
Enhanced ClinVar SNPs source for deep intronic variants with proper exon-based filtering.

This module processes ClinVar VCF files to extract deep intronic variants
that are within gene boundaries but outside of exon boundaries by a specified distance.
"""

import logging
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ..core.cache_manager import CacheManager
from ..core.ensembl_client import EnsemblClient
from .clinvar_snps import (
    _convert_to_snp_format,
    _extract_gene_symbol,
    _filter_pathogenic_variants,
    _filter_variants_by_gene_panel,
    _get_clinvar_vcf_file,
    _parse_clinvar_vcf,
    _standardize_chromosome,
)

logger = logging.getLogger(__name__)


def fetch_clinvar_snps(
    config: dict[str, Any],
    gene_panel: Optional[pd.DataFrame] = None,
    ensembl_client: Optional[EnsemblClient] = None,
) -> pd.DataFrame:
    """
    Fetch deep intronic ClinVar SNPs for the final gene panel.

    Args:
        config: Configuration dictionary containing ClinVar settings
        gene_panel: DataFrame with final gene panel (must have coordinates)
        ensembl_client: Optional Ensembl client for fetching exon data

    Returns:
        DataFrame with ClinVar SNP data filtered to deep intronic variants only
    """
    logger.info(
        "Starting ClinVar SNPs processing with proper deep intronic filtering..."
    )

    clinvar_config = config.get("snp_processing", {}).get("deep_intronic_clinvar", {})

    if not clinvar_config.get("enabled", False):
        logger.info("ClinVar SNPs processing is disabled")
        return pd.DataFrame()

    # Check if we have gene panel data for filtering
    if gene_panel is None or gene_panel.empty:
        logger.warning(
            "No gene panel provided for ClinVar filtering - skipping ClinVar processing"
        )
        return pd.DataFrame()

    # Filter to only genes included in the final panel
    if "include_in_panel" in gene_panel.columns:
        panel_genes = gene_panel[gene_panel["include_in_panel"]].copy()
    else:
        # If no include_in_panel column, assume all genes are included
        panel_genes = gene_panel.copy()

    if panel_genes.empty:
        logger.warning("No genes included in panel - skipping ClinVar processing")
        return pd.DataFrame()

    # Filter out genes with missing coordinate data
    genes_with_coords = panel_genes.dropna(
        subset=["chromosome", "gene_start", "gene_end"]
    )
    missing_coords_count = len(panel_genes) - len(genes_with_coords)

    if missing_coords_count > 0:
        logger.info(
            f"Excluding {missing_coords_count} genes with missing coordinate data from ClinVar filtering"
        )

    if genes_with_coords.empty:
        logger.warning(
            "No genes with coordinate data available - skipping ClinVar processing"
        )
        return pd.DataFrame()

    panel_genes = genes_with_coords

    logger.info(f"Processing ClinVar variants for {len(panel_genes)} panel genes")

    # Initialize Ensembl client if not provided
    if ensembl_client is None:
        cache_manager = CacheManager(
            cache_dir=str(Path(clinvar_config.get("cache_dir", ".cache/clinvar"))),
        )
        ensembl_client = EnsemblClient(cache_manager=cache_manager)

    # Get target chromosomes from panel genes to focus VCF processing
    target_chromosomes = set()
    for chrom in panel_genes["chromosome"].unique():
        if pd.notna(chrom):  # Skip NaN values
            standardized_chrom = _standardize_chromosome(str(chrom))
            if (
                standardized_chrom and standardized_chrom != "None"
            ):  # Skip None/invalid chromosomes
                target_chromosomes.add(standardized_chrom)

    logger.info(
        f"Target chromosomes for ClinVar processing: {sorted(target_chromosomes)}"
    )

    # Get or download ClinVar VCF file
    vcf_file = _get_clinvar_vcf_file(clinvar_config)
    if not vcf_file or not vcf_file.exists():
        logger.error("Could not obtain ClinVar VCF file")
        return pd.DataFrame()

    try:
        # Parse ClinVar VCF file focusing on target chromosomes
        clinvar_variants = _parse_clinvar_vcf(vcf_file, target_chromosomes)

        if clinvar_variants.empty:
            logger.warning("No variants found in ClinVar VCF file")
            return pd.DataFrame()

        logger.info(f"Parsed {len(clinvar_variants)} variants from ClinVar VCF")

        # Filter for pathogenic/likely pathogenic variants
        pathogenic_variants = _filter_pathogenic_variants(clinvar_variants)

        if pathogenic_variants.empty:
            logger.warning("No pathogenic variants found in ClinVar")
            return pd.DataFrame()

        logger.info(f"Found {len(pathogenic_variants)} pathogenic ClinVar variants")

        # Filter variants to only those within gene panel regions
        panel_filtered_variants = _filter_variants_by_gene_panel(
            pathogenic_variants, panel_genes, clinvar_config
        )

        if panel_filtered_variants.empty:
            logger.warning("No pathogenic variants found within gene panel regions")
            return pd.DataFrame()

        logger.info(
            f"Found {len(panel_filtered_variants)} pathogenic variants within gene panel regions"
        )

        # Fetch exon coordinates for panel genes
        logger.info("Fetching exon coordinates for deep intronic filtering...")
        gene_exons = _fetch_gene_exons(panel_genes, ensembl_client)

        # Apply deep intronic filtering using exon coordinates
        deep_intronic_variants = _filter_deep_intronic_with_exons(
            panel_filtered_variants, panel_genes, gene_exons, clinvar_config
        )

        if deep_intronic_variants.empty:
            logger.warning("No deep intronic variants found in gene panel regions")
            return pd.DataFrame()

        logger.info(
            f"Found {len(deep_intronic_variants)} deep intronic variants in gene panel"
        )

        # Convert to standardized SNP format
        snp_data = _convert_to_snp_format(deep_intronic_variants)

        logger.info(
            f"Successfully processed {len(snp_data)} deep intronic ClinVar SNPs"
        )
        return snp_data

    except Exception as e:
        logger.error(f"Error processing ClinVar VCF: {e}")
        import traceback

        traceback.print_exc()
        return pd.DataFrame()


def _fetch_gene_exons(
    panel_genes: pd.DataFrame, ensembl_client: EnsemblClient
) -> dict[str, list[dict[str, Any]]]:
    """
    Fetch exon coordinates for all genes in the panel.

    Args:
        panel_genes: DataFrame with gene panel
        ensembl_client: Ensembl client for API calls

    Returns:
        Dictionary mapping gene symbol to list of exon coordinates
    """
    gene_exons = {}

    # Process genes with transcript IDs first (faster)
    genes_with_transcripts = panel_genes[
        panel_genes["canonical_transcript"].notna()
    ].copy()
    genes_without_transcripts = panel_genes[
        panel_genes["canonical_transcript"].isna()
    ].copy()

    logger.info(
        f"Fetching exons for {len(genes_with_transcripts)} genes with transcript IDs..."
    )

    # Batch process genes with transcript IDs
    for _, gene in genes_with_transcripts.iterrows():
        symbol = gene["approved_symbol"]
        transcript_id = gene["canonical_transcript"]

        if pd.notna(transcript_id):
            exons = ensembl_client.get_transcript_exons(transcript_id)
            if exons:
                gene_exons[symbol] = exons
                logger.debug(
                    f"Fetched {len(exons)} exons for {symbol} from transcript {transcript_id}"
                )

    # For genes without transcript IDs, try to fetch by gene ID
    if not genes_without_transcripts.empty:
        logger.info(
            f"Fetching exons for {len(genes_without_transcripts)} genes without transcript IDs..."
        )

        for _, gene in genes_without_transcripts.iterrows():
            symbol = gene["approved_symbol"]
            gene_id = gene.get("gene_id")

            if pd.notna(gene_id):
                # Try to get canonical transcript for the gene
                gene_data = ensembl_client._make_request(
                    f"lookup/id/{gene_id}?expand=1"
                )
                if gene_data and isinstance(gene_data, dict):
                    exons = ensembl_client.get_gene_exons_by_transcript_type(
                        gene_data, "canonical"
                    )
                    if exons:
                        gene_exons[symbol] = exons
                        logger.debug(
                            f"Fetched {len(exons)} exons for {symbol} from gene {gene_id}"
                        )

    logger.info(f"Successfully fetched exon data for {len(gene_exons)} genes")
    return gene_exons


def _filter_deep_intronic_with_exons(
    variants_df: pd.DataFrame,
    panel_genes: pd.DataFrame,
    gene_exons: dict[str, list[dict[str, Any]]],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Filter variants to only deep intronic ones using actual exon coordinates.

    Deep intronic variants are defined as variants that are:
    1. Within gene boundaries (intragenic)
    2. At least 'intronic_padding' base pairs away from any exon

    Args:
        variants_df: DataFrame with ClinVar variants within gene panel regions
        panel_genes: DataFrame with final gene panel
        gene_exons: Dictionary mapping gene symbols to exon coordinates
        config: ClinVar configuration dictionary

    Returns:
        DataFrame with only deep intronic variants
    """
    if variants_df.empty:
        return pd.DataFrame()

    # Get intronic padding from config (default 50bp)
    intronic_padding = config.get("intronic_padding", 50)
    logger.info(
        f"Applying deep intronic filter with {intronic_padding}bp padding from exon boundaries"
    )

    # Prepare gene coordinate lookup for efficiency
    gene_coords: dict[str, list[dict[str, Any]]] = {}
    for _, gene in panel_genes.iterrows():
        gene_chr = _standardize_chromosome(str(gene["chromosome"]))
        gene_symbol = gene["approved_symbol"]

        if pd.isna(gene["gene_start"]) or pd.isna(gene["gene_end"]):
            continue

        if gene_chr not in gene_coords:
            gene_coords[gene_chr] = []

        gene_coords[gene_chr].append(
            {
                "symbol": gene_symbol,
                "start": int(gene["gene_start"]),
                "end": int(gene["gene_end"]),
            }
        )

    # Pre-process exon data for efficient lookup
    gene_exon_ranges = {}
    for gene_symbol, exons in gene_exons.items():
        if not exons:
            continue

        # Create list of exon ranges with padding
        exon_ranges = []
        for exon in exons:
            if exon.get("start") is not None and exon.get("end") is not None:
                exon_ranges.append(
                    {
                        "start": exon["start"] - intronic_padding,
                        "end": exon["end"] + intronic_padding,
                    }
                )

        if exon_ranges:
            # Sort by start position for efficient searching
            exon_ranges.sort(key=lambda x: x["start"])
            gene_exon_ranges[gene_symbol] = exon_ranges

    # Process variants in batches by chromosome
    deep_intronic_mask = []
    total_variants = len(variants_df)
    processed = 0

    for _, variant in variants_df.iterrows():
        if processed % 1000 == 0 and processed > 0:
            logger.debug(
                f"Processed {processed}/{total_variants} variants for deep intronic filtering"
            )
        processed += 1

        variant_chr = _standardize_chromosome(str(variant["chromosome"]))
        variant_pos = variant["position"]
        variant_gene = _extract_gene_symbol(variant.get("geneinfo", ""))

        # Quick check: if no genes on this chromosome, skip
        if variant_chr not in gene_coords:
            deep_intronic_mask.append(False)
            continue

        # Find overlapping genes
        overlapping_genes = []
        for gene_info in gene_coords[variant_chr]:
            if gene_info["start"] <= variant_pos <= gene_info["end"]:
                # If variant has gene annotation, check if it matches
                if variant_gene and variant_gene != gene_info["symbol"]:
                    continue
                overlapping_genes.append(gene_info["symbol"])

        if not overlapping_genes:
            deep_intronic_mask.append(False)
            continue

        # Check if variant is deep intronic for any overlapping gene
        is_deep_intronic = False

        for gene_symbol in overlapping_genes:
            if gene_symbol not in gene_exon_ranges:
                # No exon data for this gene - skip it
                continue

            exon_ranges = gene_exon_ranges[gene_symbol]

            # Check if variant is within any padded exon range
            is_near_exon = False
            for exon_range in exon_ranges:
                if exon_range["start"] <= variant_pos <= exon_range["end"]:
                    is_near_exon = True
                    break
                # Since exons are sorted, can break early if past variant position
                if exon_range["start"] > variant_pos:
                    break

            # If not near any exon (including padding), it's deep intronic
            if not is_near_exon:
                is_deep_intronic = True
                break

        deep_intronic_mask.append(is_deep_intronic)

    # Apply mask to filter variants
    result_df = variants_df[deep_intronic_mask].copy()

    logger.info(
        f"Identified {len(result_df)} truly deep intronic variants (>{intronic_padding}bp from exons)"
    )

    return result_df
