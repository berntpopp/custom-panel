"""
Highly optimized ClinVar SNPs source using vectorized operations and parallel processing.

This module processes ClinVar VCF files to extract deep intronic variants
using efficient pandas operations and parallel processing.
"""

import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..core.cache_manager import CacheManager
from ..core.ensembl_client import EnsemblClient

logger = logging.getLogger(__name__)


def fetch_clinvar_snps(config: dict, gene_panel: Optional[pd.DataFrame] = None,
                      ensembl_client: Optional[EnsemblClient] = None) -> pd.DataFrame:
    """
    Fetch deep intronic ClinVar SNPs for the final gene panel using vectorized operations.

    Args:
        config: Configuration dictionary containing ClinVar settings
        gene_panel: DataFrame with final gene panel (must have coordinates)
        ensembl_client: Optional Ensembl client for fetching exon data

    Returns:
        DataFrame with ClinVar SNP data filtered to deep intronic variants only
    """
    logger.info("Starting optimized ClinVar SNPs processing with vectorized deep intronic filtering...")

    clinvar_config = config.get("snp_processing", {}).get("deep_intronic_clinvar", {})

    if not clinvar_config.get("enabled", False):
        logger.info("ClinVar SNPs processing is disabled")
        return pd.DataFrame()

    # Check if we have gene panel data for filtering
    if gene_panel is None or gene_panel.empty:
        logger.warning("No gene panel provided for ClinVar filtering - skipping ClinVar processing")
        return pd.DataFrame()

    # Filter to only genes included in the final panel
    if 'include_in_panel' in gene_panel.columns:
        panel_genes = gene_panel[gene_panel['include_in_panel'] == True].copy()
    else:
        # If no include_in_panel column, assume all genes are included
        panel_genes = gene_panel.copy()

    if panel_genes.empty:
        logger.warning("No genes included in panel - skipping ClinVar processing")
        return pd.DataFrame()

    # Filter out genes with missing coordinate data
    genes_with_coords = panel_genes.dropna(subset=['chromosome', 'gene_start', 'gene_end'])
    missing_coords_count = len(panel_genes) - len(genes_with_coords)

    if missing_coords_count > 0:
        logger.info(f"Excluding {missing_coords_count} genes with missing coordinate data from ClinVar filtering")

    if genes_with_coords.empty:
        logger.warning("No genes with coordinate data available - skipping ClinVar processing")
        return pd.DataFrame()

    panel_genes = genes_with_coords

    logger.info(f"Processing ClinVar variants for {len(panel_genes)} panel genes")

    # Initialize Ensembl client if not provided
    if ensembl_client is None:
        cache_manager = CacheManager(
            cache_dir=Path(clinvar_config.get("cache_dir", ".cache/clinvar")),
            default_ttl=clinvar_config.get("cache_expiry_days", 30) * 24 * 60 * 60
        )
        ensembl_client = EnsemblClient(cache_manager=cache_manager)

    # Get target chromosomes from panel genes to focus VCF processing
    target_chromosomes = set()
    for chrom in panel_genes['chromosome'].unique():
        if pd.notna(chrom):  # Skip NaN values
            standardized_chrom = _standardize_chromosome(str(chrom))
            if standardized_chrom and standardized_chrom != 'None':  # Skip None/invalid chromosomes
                target_chromosomes.add(standardized_chrom)

    logger.info(f"Target chromosomes for ClinVar processing: {sorted(target_chromosomes)}")

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

        # Filter variants to only those within gene panel regions (vectorized)
        panel_filtered_variants = _filter_variants_by_gene_panel_vectorized(
            pathogenic_variants, panel_genes, clinvar_config
        )

        if panel_filtered_variants.empty:
            logger.warning("No pathogenic variants found within gene panel regions")
            return pd.DataFrame()

        logger.info(f"Found {len(panel_filtered_variants)} pathogenic variants within gene panel regions")

        # Fetch exon coordinates for panel genes
        logger.info("Fetching exon coordinates for deep intronic filtering...")
        gene_exons = _fetch_gene_exons(panel_genes, ensembl_client)

        # Apply deep intronic filtering using vectorized operations
        deep_intronic_variants = _filter_deep_intronic_vectorized(
            panel_filtered_variants, panel_genes, gene_exons, clinvar_config
        )

        if deep_intronic_variants.empty:
            logger.warning("No deep intronic variants found in gene panel regions")
            return pd.DataFrame()

        logger.info(f"Found {len(deep_intronic_variants)} deep intronic variants in gene panel")

        # Convert to standardized SNP format
        snp_data = _convert_to_snp_format(deep_intronic_variants)

        logger.info(f"Successfully processed {len(snp_data)} deep intronic ClinVar SNPs")
        return snp_data

    except Exception as e:
        logger.error(f"Error processing ClinVar VCF: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()


def _filter_variants_by_gene_panel_vectorized(
    variants_df: pd.DataFrame,
    panel_genes: pd.DataFrame,
    config: dict
) -> pd.DataFrame:
    """
    Filter ClinVar variants to only those within gene panel genomic regions using vectorized operations.
    """
    if variants_df.empty or panel_genes.empty:
        return pd.DataFrame()

    logger.info(f"Filtering {len(variants_df)} variants against {len(panel_genes)} genes using vectorized operations")

    # Standardize chromosomes
    variants_df['std_chr'] = variants_df['chromosome'].apply(lambda x: _standardize_chromosome(str(x)))
    variants_df['gene_symbol'] = variants_df['geneinfo'].apply(_extract_gene_symbol)

    panel_genes['std_chr'] = panel_genes['chromosome'].apply(lambda x: _standardize_chromosome(str(x)))

    # Convert coordinates to numeric
    panel_genes['gene_start'] = pd.to_numeric(panel_genes['gene_start'], errors='coerce')
    panel_genes['gene_end'] = pd.to_numeric(panel_genes['gene_end'], errors='coerce')
    variants_df['position'] = pd.to_numeric(variants_df['position'], errors='coerce')

    # Process by chromosome for efficiency
    filtered_dfs = []

    for chrom in panel_genes['std_chr'].unique():
        if pd.isna(chrom) or chrom == 'None':
            continue

        # Get data for this chromosome
        chrom_variants = variants_df[variants_df['std_chr'] == chrom].copy()
        chrom_genes = panel_genes[panel_genes['std_chr'] == chrom].copy()

        if chrom_variants.empty or chrom_genes.empty:
            continue

        # For each variant, check if it falls within any gene using vectorized operations
        # Create a cross join effect by repeating variants for each gene
        n_variants = len(chrom_variants)
        n_genes = len(chrom_genes)

        # Repeat variant positions for comparison with all genes
        variant_positions = np.repeat(chrom_variants['position'].values, n_genes)
        variant_indices = np.repeat(chrom_variants.index.values, n_genes)

        # Tile gene boundaries to match variant repetitions
        gene_starts = np.tile(chrom_genes['gene_start'].values, n_variants)
        gene_ends = np.tile(chrom_genes['gene_end'].values, n_variants)

        # Vectorized comparison
        within_gene = (variant_positions >= gene_starts) & (variant_positions <= gene_ends)

        # Get unique variant indices that fall within at least one gene
        valid_variant_indices = np.unique(variant_indices[within_gene])

        if len(valid_variant_indices) > 0:
            filtered_dfs.append(chrom_variants.loc[valid_variant_indices])

    if filtered_dfs:
        result_df = pd.concat(filtered_dfs, ignore_index=True)
        # Remove temporary columns
        result_df = result_df.drop(columns=['std_chr', 'gene_symbol'], errors='ignore')
        logger.info(f"Filtered to {len(result_df)} variants within gene panel regions")
        return result_df
    else:
        return pd.DataFrame()


def _filter_deep_intronic_vectorized(
    variants_df: pd.DataFrame,
    panel_genes: pd.DataFrame,
    gene_exons: dict[str, list[dict[str, Any]]],
    config: dict
) -> pd.DataFrame:
    """
    Filter variants to only deep intronic ones using vectorized operations.
    """
    if variants_df.empty:
        return pd.DataFrame()

    # Get intronic padding from config (default 50bp)
    intronic_padding = config.get('intronic_padding', 50)
    logger.info(f"Applying vectorized deep intronic filter with {intronic_padding}bp padding")

    # Prepare exon intervals for all genes
    all_exon_intervals = []
    for gene_symbol, exons in gene_exons.items():
        for exon in exons:
            if exon.get('start') is not None and exon.get('end') is not None:
                all_exon_intervals.append({
                    'gene_symbol': gene_symbol,
                    'chromosome': _standardize_chromosome(str(exon.get('chromosome', ''))),
                    'exon_start': exon['start'] - intronic_padding,
                    'exon_end': exon['end'] + intronic_padding
                })

    if not all_exon_intervals:
        logger.warning("No valid exon intervals found")
        return pd.DataFrame()

    exon_df = pd.DataFrame(all_exon_intervals)

    # Standardize variant chromosomes
    variants_df['std_chr'] = variants_df['chromosome'].apply(lambda x: _standardize_chromosome(str(x)))
    variants_df['gene_symbol'] = variants_df['geneinfo'].apply(_extract_gene_symbol)

    # Use parallel processing for chromosome batches
    n_cores = min(mp.cpu_count(), 8)  # Limit to 8 cores
    chromosomes = variants_df['std_chr'].unique()

    # Split chromosomes into batches for parallel processing
    chrom_batches = np.array_split(chromosomes, n_cores)

    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        futures = []

        for batch in chrom_batches:
            if len(batch) == 0:
                continue

            batch_variants = variants_df[variants_df['std_chr'].isin(batch)].copy()
            batch_exons = exon_df[exon_df['chromosome'].isin(batch)].copy()
            batch_genes = panel_genes[panel_genes['chromosome'].apply(
                lambda x: _standardize_chromosome(str(x))
            ).isin(batch)].copy()

            future = executor.submit(
                _process_chromosome_batch,
                batch_variants,
                batch_exons,
                batch_genes,
                intronic_padding
            )
            futures.append(future)

        # Collect results
        deep_intronic_indices = []
        for future in as_completed(futures):
            try:
                indices = future.result()
                deep_intronic_indices.extend(indices)
            except Exception as e:
                logger.error(f"Error in parallel processing: {e}")

    if deep_intronic_indices:
        result_df = variants_df.loc[deep_intronic_indices].copy()
        result_df = result_df.drop(columns=['std_chr', 'gene_symbol'], errors='ignore')
        logger.info(f"Identified {len(result_df)} truly deep intronic variants")
        return result_df
    else:
        return pd.DataFrame()


def _process_chromosome_batch(
    variants_df: pd.DataFrame,
    exon_df: pd.DataFrame,
    panel_genes: pd.DataFrame,
    intronic_padding: int
) -> list[int]:
    """
    Process a batch of chromosomes to identify deep intronic variants.
    """
    deep_intronic_indices = []

    for idx, variant in variants_df.iterrows():
        variant_pos = variant['position']
        variant_chr = variant['std_chr']
        variant_gene = variant['gene_symbol']

        # Find overlapping genes
        gene_matches = panel_genes[
            (panel_genes['chromosome'].apply(lambda x: _standardize_chromosome(str(x))) == variant_chr) &
            (panel_genes['gene_start'] <= variant_pos) &
            (panel_genes['gene_end'] >= variant_pos)
        ]

        if gene_matches.empty:
            continue

        # Check each overlapping gene
        is_deep_intronic = False
        for _, gene in gene_matches.iterrows():
            gene_symbol = gene['approved_symbol']

            # Get exons for this gene
            gene_exons = exon_df[
                (exon_df['gene_symbol'] == gene_symbol) &
                (exon_df['chromosome'] == variant_chr)
            ]

            if gene_exons.empty:
                continue

            # Check if variant is outside all padded exon boundaries
            is_outside_all_exons = True
            for _, exon in gene_exons.iterrows():
                if exon['exon_start'] <= variant_pos <= exon['exon_end']:
                    is_outside_all_exons = False
                    break

            if is_outside_all_exons:
                is_deep_intronic = True
                break

        if is_deep_intronic:
            deep_intronic_indices.append(idx)

    return deep_intronic_indices


# Import helper functions from enhanced implementation
# Import remaining functions from original implementation
from .clinvar_snps import (
    _convert_to_snp_format,
    _extract_gene_symbol,
    _filter_pathogenic_variants,
    _get_clinvar_vcf_file,
    _parse_clinvar_vcf,
    _standardize_chromosome,
)
from .clinvar_snps_enhanced import _fetch_gene_exons
