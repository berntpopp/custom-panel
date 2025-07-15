"""
Highly optimized ClinVar SNPs source using vectorized operations and parallel processing.

This module processes ClinVar VCF files to extract deep intronic variants
using efficient pandas operations and parallel processing.
"""

from __future__ import annotations

import gzip
import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
import requests

from ..core.cache_manager import CacheManager
from ..core.ensembl_client import EnsemblClient

logger = logging.getLogger(__name__)


def fetch_clinvar_snps(
    config: dict[str, Any],
    gene_panel: Optional[pd.DataFrame] = None,
    ensembl_client: Optional[EnsemblClient] = None,
    harmonizer: Optional[Any] = None,
) -> pd.DataFrame:
    """
    Fetch deep intronic ClinVar SNPs for the final gene panel using vectorized operations.

    Args:
        config: Configuration dictionary containing ClinVar settings
        gene_panel: DataFrame with final gene panel (must have coordinates)
        ensembl_client: Optional Ensembl client for fetching exon data

    Returns:
        DataFrame with ClinVar SNP data filtered to deep intronic variants only
    """
    logger.info(
        "Starting optimized ClinVar SNPs processing with vectorized deep intronic filtering..."
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
        panel_genes = gene_panel.copy()

    if panel_genes.empty:
        logger.warning("No genes marked for inclusion in panel")
        return pd.DataFrame()

    logger.info(f"Processing ClinVar SNPs for {len(panel_genes)} genes")

    # Check if panel genes have required coordinate columns
    required_coords = ["chromosome", "gene_start", "gene_end"]
    if not all(col in panel_genes.columns for col in required_coords):
        logger.warning(
            f"Gene panel missing required coordinate columns {required_coords}"
        )
        logger.warning(f"Available columns: {list(panel_genes.columns)}")
        return pd.DataFrame()

    # Remove genes without coordinate data
    panel_genes = panel_genes.dropna(subset=required_coords)
    if panel_genes.empty:
        logger.warning("No genes with valid coordinate data")
        return pd.DataFrame()

    # Create Ensembl client if not provided
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

        # Filter out large CNVs and structural variants (keep only small indels up to 50bp)
        size_filtered_variants = _filter_variant_size(
            pathogenic_variants, clinvar_config
        )

        if size_filtered_variants.empty:
            logger.warning("No variants found after size filtering")
            return pd.DataFrame()

        logger.info(
            f"Found {len(size_filtered_variants)} variants after size filtering (removed large CNVs)"
        )

        # Filter variants to only those within gene panel regions using vectorized operations
        panel_filtered_variants = _filter_variants_by_gene_panel_vectorized(
            size_filtered_variants, panel_genes, clinvar_config
        )

        if panel_filtered_variants.empty:
            logger.warning("No pathogenic variants found within gene panel regions")
            return pd.DataFrame()

        logger.info(
            f"Found {len(panel_filtered_variants)} pathogenic variants within gene panel regions"
        )

        # Apply deep intronic filtering using exon coordinates
        deep_intronic_variants = _filter_deep_intronic_vectorized(
            panel_filtered_variants, panel_genes, ensembl_client, clinvar_config
        )

        if deep_intronic_variants.empty:
            logger.warning("No deep intronic variants found after filtering")
            return pd.DataFrame()

        logger.info(
            f"Found {len(deep_intronic_variants)} deep intronic ClinVar variants"
        )

        # Convert to standardized SNP format
        snp_formatted = _convert_to_snp_format(deep_intronic_variants, harmonizer)

        logger.info(f"Successfully processed {len(snp_formatted)} ClinVar SNPs")
        return snp_formatted

    except Exception as e:
        logger.error(f"Error processing ClinVar SNPs: {str(e)}")
        return pd.DataFrame()


def _filter_variants_by_gene_panel_vectorized(
    variants_df: pd.DataFrame, panel_genes: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """
    Filter ClinVar variants to only those within gene panel genomic regions using vectorized operations.

    Args:
        variants_df: DataFrame with ClinVar variants
        panel_genes: DataFrame with final gene panel (must have genomic coordinates)
        config: ClinVar configuration dictionary

    Returns:
        Filtered DataFrame with variants only within gene panel regions
    """
    if variants_df.empty or panel_genes.empty:
        return pd.DataFrame()

    logger.info(
        f"Filtering {len(variants_df)} variants against {len(panel_genes)} genes using vectorized operations"
    )

    # Prepare variants with standardized chromosomes and gene symbols
    variants_prepared = variants_df.copy()
    variants_prepared["std_chr"] = variants_prepared["chromosome"].apply(
        lambda x: _standardize_chromosome(str(x))
    )
    variants_prepared["gene_symbol"] = variants_prepared["geneinfo"].apply(
        _extract_gene_symbol
    )

    # Prepare genes with standardized chromosomes
    genes_prepared = panel_genes.copy()
    genes_prepared["std_chr"] = genes_prepared["chromosome"].apply(
        lambda x: _standardize_chromosome(str(x))
    )

    # Use vectorized operations for chromosome-position filtering
    filtered_variants = []

    # Process by chromosome for efficiency
    for chrom in genes_prepared["std_chr"].unique():
        if pd.isna(chrom) or chrom == "None":
            continue

        # Get variants and genes for this chromosome
        chrom_variants = variants_prepared[variants_prepared["std_chr"] == chrom]
        chrom_genes = genes_prepared[genes_prepared["std_chr"] == chrom]

        if chrom_variants.empty or chrom_genes.empty:
            continue

        # Vectorized position filtering using NumPy broadcasting
        variant_positions = chrom_variants["position"].values
        gene_starts = chrom_genes["gene_start"].values
        gene_ends = chrom_genes["gene_end"].values

        # Create boolean mask for variants within any gene region
        # This uses broadcasting to compare all variants against all genes
        within_gene_mask = (
            (variant_positions[:, np.newaxis] >= gene_starts[np.newaxis, :])
            & (variant_positions[:, np.newaxis] <= gene_ends[np.newaxis, :])
        ).any(axis=1)

        # Filter variants using the mask
        filtered_chrom_variants = chrom_variants[within_gene_mask]
        filtered_variants.append(filtered_chrom_variants)

    # Combine results from all chromosomes
    if filtered_variants:
        result_df = pd.concat(filtered_variants, ignore_index=True)
    else:
        result_df = pd.DataFrame()

    logger.info(
        f"Found {len(result_df)} variants within gene panel regions after vectorized filtering"
    )

    return result_df


def _filter_deep_intronic_vectorized(
    variants_df: pd.DataFrame,
    panel_genes: pd.DataFrame,
    ensembl_client: EnsemblClient,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Filter variants to only deep intronic ones using vectorized operations and parallel processing.

    Args:
        variants_df: DataFrame with ClinVar variants within gene panel regions
        panel_genes: DataFrame with final gene panel
        ensembl_client: Ensembl client for fetching exon data
        config: ClinVar configuration dictionary

    Returns:
        DataFrame with only deep intronic variants
    """
    if variants_df.empty:
        return pd.DataFrame()

    # Get intronic padding from config (default 50bp)
    intronic_padding = config.get("intronic_padding", 50)
    logger.info(
        f"Applying vectorized deep intronic filter with {intronic_padding}bp padding from exon boundaries"
    )

    # Fetch exon coordinates for all genes
    gene_exons = _fetch_gene_exons(panel_genes, ensembl_client)

    if not gene_exons:
        logger.warning("No exon data available - cannot apply deep intronic filtering")
        return pd.DataFrame()

    # Group variants by chromosome for parallel processing
    chromosomes = variants_df["chromosome"].apply(_standardize_chromosome).unique()
    chromosome_groups = []

    for chrom in chromosomes:
        if pd.isna(chrom) or chrom == "None":
            continue

        chrom_variants = variants_df[
            variants_df["chromosome"].apply(_standardize_chromosome) == chrom
        ]
        if not chrom_variants.empty:
            chromosome_groups.append((chrom, chrom_variants))

    # Process chromosomes in parallel
    max_workers = min(mp.cpu_count(), 8)  # Limit to 8 cores max
    deep_intronic_variants = []

    if len(chromosome_groups) > 1 and max_workers > 1:
        logger.info(f"Processing {len(chromosome_groups)} chromosomes in parallel")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit tasks for each chromosome
            future_to_chrom = {
                executor.submit(
                    _process_chromosome_batch,
                    chrom,
                    chrom_variants,
                    gene_exons,
                    intronic_padding,
                ): chrom
                for chrom, chrom_variants in chromosome_groups
            }

            # Collect results
            for future in as_completed(future_to_chrom):
                chrom = future_to_chrom[future]
                try:
                    chrom_deep_intronic = future.result()
                    if not chrom_deep_intronic.empty:
                        deep_intronic_variants.append(chrom_deep_intronic)
                except Exception as e:
                    logger.error(f"Error processing chromosome {chrom}: {e}")
    else:
        # Sequential processing for single chromosome or limited workers
        for chrom, chrom_variants in chromosome_groups:
            chrom_deep_intronic = _process_chromosome_batch(
                chrom, chrom_variants, gene_exons, intronic_padding
            )
            if not chrom_deep_intronic.empty:
                deep_intronic_variants.append(chrom_deep_intronic)

    # Combine results
    if deep_intronic_variants:
        result_df = pd.concat(deep_intronic_variants, ignore_index=True)
    else:
        result_df = pd.DataFrame()

    logger.info(
        f"Identified {len(result_df)} deep intronic variants after vectorized filtering"
    )

    return result_df


def _process_chromosome_batch(
    chrom: str,
    chrom_variants: pd.DataFrame,
    gene_exons: dict[str, list[dict[str, Any]]],
    intronic_padding: int,
) -> pd.DataFrame:
    """
    Process a batch of variants from a single chromosome.

    Args:
        chrom: Chromosome identifier
        chrom_variants: Variants from this chromosome
        gene_exons: Dictionary mapping gene symbols to exon coordinates
        intronic_padding: Distance from exon boundaries in bp

    Returns:
        DataFrame with deep intronic variants from this chromosome
    """
    deep_intronic_variants = []

    for _, variant in chrom_variants.iterrows():
        variant_pos = variant["position"]
        variant_gene = _extract_gene_symbol(variant.get("geneinfo", ""))

        # Check if variant is deep intronic
        if variant_gene and variant_gene in gene_exons:
            exons = gene_exons[variant_gene]
            is_deep_intronic = True

            # Check distance from all exons using vectorized operations
            if exons:
                exon_starts = np.array([exon["start"] for exon in exons])
                exon_ends = np.array([exon["end"] for exon in exons])

                # Calculate minimum distance to any exon boundary
                dist_to_starts = np.abs(variant_pos - exon_starts)
                dist_to_ends = np.abs(variant_pos - exon_ends)

                # Check if variant is within any exon or too close to boundaries
                within_exon = (variant_pos >= exon_starts) & (variant_pos <= exon_ends)
                too_close_to_boundary = (dist_to_starts < intronic_padding) | (
                    dist_to_ends < intronic_padding
                )

                if within_exon.any() or too_close_to_boundary.any():
                    is_deep_intronic = False

            if is_deep_intronic:
                deep_intronic_variants.append(variant)
        else:
            # If no exon data available, keep the variant (conservative approach)
            deep_intronic_variants.append(variant)

    if deep_intronic_variants:
        return pd.DataFrame(deep_intronic_variants)
    else:
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


# Helper functions from original implementation
def _get_clinvar_vcf_file(config: dict[str, Any]) -> Optional[Path]:
    """
    Get ClinVar VCF file - either from local path or by downloading.

    Args:
        config: ClinVar configuration dictionary

    Returns:
        Path to ClinVar VCF file or None if failed
    """
    # Check if local file path is provided
    vcf_path = config.get("vcf_path")
    if vcf_path:
        vcf_file = Path(vcf_path)
        if vcf_file.exists():
            logger.info(f"Using provided ClinVar VCF file: {vcf_path}")
            return vcf_file
        else:
            logger.warning(f"Provided ClinVar VCF file not found: {vcf_path}")

    # Download ClinVar VCF if not provided or not found
    vcf_url = config.get("vcf_url")
    if not vcf_url:
        logger.error("No ClinVar VCF URL provided for download")
        return None

    cache_dir = Path(config.get("cache_dir", ".cache/clinvar"))
    cache_expiry_days = config.get("cache_expiry_days", 30)

    # Create cache directory
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Download VCF file
    return _download_clinvar_vcf(vcf_url, cache_dir, cache_expiry_days)


def _download_clinvar_vcf(
    url: str, cache_dir: Path, cache_expiry_days: int
) -> Optional[Path]:
    """
    Download ClinVar VCF file with caching.

    Args:
        url: URL to download from
        cache_dir: Directory to cache the file
        cache_expiry_days: Number of days before re-downloading

    Returns:
        Path to downloaded file or None if failed
    """
    # Determine filename from URL
    filename = url.split("/")[-1]
    if not filename.endswith((".vcf", ".vcf.gz")):
        filename = "clinvar.vcf.gz"

    cache_file = cache_dir / filename

    # Check if cached file exists and is recent enough
    if cache_file.exists():
        file_age = datetime.now() - datetime.fromtimestamp(cache_file.stat().st_mtime)
        if file_age.days < cache_expiry_days:
            logger.info(f"Using cached ClinVar VCF: {cache_file}")
            return cache_file
        else:
            logger.info(
                f"Cached ClinVar VCF is {file_age.days} days old, re-downloading..."
            )

    # Download the file
    logger.info(f"Downloading ClinVar VCF from: {url}")
    logger.info("This may take several minutes due to file size (~500MB)")

    try:
        # Stream download for large files
        response = requests.get(url, stream=True, timeout=300)
        response.raise_for_status()

        # Get file size for progress reporting
        total_size = int(response.headers.get("content-length", 0))
        downloaded_size = 0

        with open(cache_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded_size += len(chunk)

                    # Progress reporting every 10MB
                    if total_size > 0 and downloaded_size % (10 * 1024 * 1024) == 0:
                        progress = (downloaded_size / total_size) * 100
                        logger.info(
                            f"Downloaded {progress:.1f}% ({downloaded_size / 1024 / 1024:.1f}MB)"
                        )

        logger.info(f"Successfully downloaded ClinVar VCF to: {cache_file}")
        return cache_file

    except Exception as e:
        logger.error(f"Error downloading ClinVar VCF: {str(e)}")
        return None


def _parse_clinvar_vcf(
    vcf_file: Path, target_chromosomes: Optional[set[str]] = None
) -> pd.DataFrame:
    """
    Parse ClinVar VCF file to extract variant information.

    Args:
        vcf_file: Path to ClinVar VCF file (can be gzipped)
        target_chromosomes: Set of chromosomes to focus on (for efficiency)

    Returns:
        DataFrame with parsed variant data
    """
    variants = []

    # Determine if file is gzipped
    open_func = gzip.open if vcf_file.suffix == ".gz" else open
    mode = "rt" if vcf_file.suffix == ".gz" else "r"

    with open_func(vcf_file, mode) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("#"):
                continue  # Skip header lines

            try:
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom = fields[0]

                # Skip chromosomes not in our target set (if specified)
                if (
                    target_chromosomes
                    and _standardize_chromosome(chrom) not in target_chromosomes
                ):
                    continue

                pos = int(fields[1])
                variant_id = fields[2]  # Usually rs number
                ref = fields[3]
                alt = fields[4]
                info = fields[7]

                # Parse INFO field for ClinVar-specific information
                info_dict = _parse_info_field(info)

                variant = {
                    "chromosome": chrom,
                    "position": pos,
                    "rsid": variant_id if variant_id.startswith("rs") else None,
                    "ref_allele": ref,
                    "alt_allele": alt,
                    "clnsig": info_dict.get("CLNSIG", ""),
                    "geneinfo": info_dict.get("GENEINFO", ""),
                    "clndn": info_dict.get("CLNDN", ""),
                    "clnrevstat": info_dict.get("CLNREVSTAT", ""),
                    "clnvc": info_dict.get("CLNVC", ""),
                    "origin": info_dict.get("ORIGIN", ""),
                    "clnvcso": info_dict.get("CLNVCSO", ""),
                    "clnvi": info_dict.get("CLNVI", ""),
                }

                variants.append(variant)

                # Optional memory management for extremely large files
                # Remove artificial limit to allow processing of all variants

                # Progress reporting
                if line_num % 100000 == 0:
                    logger.info(
                        f"Processed {line_num} lines, found {len(variants)} variants"
                    )

            except (ValueError, IndexError) as e:
                logger.warning(f"Error parsing line {line_num}: {str(e)}")
                continue

    logger.info(f"Parsed {len(variants)} variants from ClinVar VCF")
    return pd.DataFrame(variants)


def _parse_info_field(info: str) -> dict[str, Any]:
    """
    Parse VCF INFO field into dictionary.

    Args:
        info: INFO field string from VCF

    Returns:
        Dictionary of INFO field values
    """
    info_dict = {}

    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = "True"

    return info_dict


def _filter_pathogenic_variants(variants_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter for pathogenic and likely pathogenic variants.

    Args:
        variants_df: DataFrame with ClinVar variants

    Returns:
        Filtered DataFrame with pathogenic variants
    """
    if variants_df.empty:
        return variants_df

    # ClinVar clinical significance values (text format)
    pathogenic_values = [
        "Pathogenic",
        "Likely_pathogenic",
        "Pathogenic/Likely_pathogenic",
        "Likely_pathogenic/Pathogenic",
    ]

    # Filter for exact pathogenic values
    pathogenic_mask = variants_df["clnsig"].isin(pathogenic_values)

    # Also include variants with pathogenic keywords in complex classifications
    pathogenic_text_mask = variants_df["clnsig"].str.contains(
        r"\bPathogenic\b|\bLikely_pathogenic\b", case=False, na=False, regex=True
    )

    # Combine filters
    combined_mask = pathogenic_mask | pathogenic_text_mask

    filtered_df = variants_df[combined_mask].copy()

    logger.info(
        f"Filtered to {len(filtered_df)} pathogenic variants from {len(variants_df)} total variants"
    )

    return filtered_df


def _filter_variant_size(
    variants_df: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """
    Filter variants by size to exclude large CNVs and structural variants.

    Args:
        variants_df: DataFrame with ClinVar variants
        config: ClinVar configuration dictionary

    Returns:
        Filtered DataFrame with variants <= max_indel_size
    """
    if variants_df.empty:
        return variants_df

    # Get maximum indel size from config (default 50bp)
    max_indel_size = config.get("max_indel_size", 50)

    logger.info(f"Filtering variants by size (max indel size: {max_indel_size}bp)")

    # Calculate variant size as the maximum of ref and alt allele lengths
    ref_lengths = variants_df["ref_allele"].astype(str).str.len()
    alt_lengths = variants_df["alt_allele"].astype(str).str.len()

    # For indels, the size is the absolute difference between ref and alt lengths
    # For SNVs, both lengths are 1, so difference is 0
    variant_sizes = np.maximum(ref_lengths, alt_lengths)

    # Filter variants by size
    size_mask = variant_sizes <= max_indel_size

    # Log some statistics about filtered variants
    large_variants_count = (~size_mask).sum()
    large_variants_examples = variants_df[~size_mask].head(3)

    if large_variants_count > 0:
        logger.info(
            f"Filtering out {large_variants_count} large variants (>{max_indel_size}bp)"
        )

        # Log examples of filtered variants
        for _, variant in large_variants_examples.iterrows():
            ref_len = len(str(variant["ref_allele"]))
            alt_len = len(str(variant["alt_allele"]))
            max_len = max(ref_len, alt_len)
            logger.debug(
                f"Filtered large variant: {variant['chromosome']}:{variant['position']} "
                f"(ref={ref_len}bp, alt={alt_len}bp, max={max_len}bp) "
                f"- {variant.get('geneinfo', 'Unknown gene')}"
            )

    filtered_df = variants_df[size_mask].copy()

    logger.info(
        f"Kept {len(filtered_df)} variants after size filtering "
        f"(removed {large_variants_count} large variants)"
    )

    return filtered_df


def _convert_to_snp_format(
    variants_df: pd.DataFrame, harmonizer: Optional[Any] = None
) -> pd.DataFrame:
    """
    Convert ClinVar variants to standardized SNP format with harmonization.

    Args:
        variants_df: DataFrame with ClinVar variants
        harmonizer: Optional SNP harmonizer for ID resolution and coordinate liftover

    Returns:
        DataFrame in standardized SNP format compatible with HTML reporting
    """
    if variants_df.empty:
        return pd.DataFrame()

    # Extract gene symbols from GENEINFO field
    variants_df["gene_symbol"] = variants_df["geneinfo"].apply(_extract_gene_symbol)

    # Create allele string format (ref/alt)
    allele_strings = (
        variants_df["ref_allele"].astype(str)
        + "/"
        + variants_df["alt_allele"].astype(str)
    )

    # Create initial SNP data with hg38 coordinates from VCF
    snp_data = pd.DataFrame(
        {
            # Standard SNP identification columns
            "snp": variants_df["rsid"],  # Use rsID from VCF when available
            "source": "ClinVar",
            "category": "deep_intronic_clinvar",
            "snp_type": "deep_intronic_clinvar",
            # hg38 coordinate columns (from VCF)
            "hg38_chromosome": variants_df["chromosome"].apply(_standardize_chromosome),
            "hg38_start": variants_df["position"],
            "hg38_end": variants_df["position"],  # SNVs have same start/end
            "hg38_strand": 1,  # Default positive strand
            "hg38_allele_string": allele_strings,
            # hg19 coordinate columns (will be populated by harmonizer)
            "hg19_chromosome": "",
            "hg19_start": "",
            "hg19_end": "",
            "hg19_strand": "",
            "hg19_allele_string": "",
            # Additional ClinVar-specific columns
            "clinical_significance": variants_df["clnsig"],
            "gene_symbol": variants_df["gene_symbol"],
            "condition": variants_df["clndn"],
            "review_status": variants_df["clnrevstat"],
            "variant_type": variants_df["clnvc"],
            "origin": variants_df["origin"],
            "molecular_consequence": variants_df["clnvcso"],
            "variant_id": variants_df["clnvi"],
            # Metadata columns
            "processing_date": datetime.now().strftime("%Y-%m-%d"),
            "source_file": "ClinVar VCF",
        }
    )

    # Apply harmonization if harmonizer is provided
    if harmonizer is not None:
        try:
            logger.info(f"Harmonizing {len(snp_data)} ClinVar SNPs")

            # Handle missing rsIDs by creating coordinate-based IDs for harmonization
            missing_rsid_mask = snp_data["snp"].isna()
            if missing_rsid_mask.any():
                # Create temporary coordinate-based IDs for harmonization
                coord_ids = (
                    snp_data.loc[missing_rsid_mask, "hg38_chromosome"].astype(str)
                    + ":"
                    + snp_data.loc[missing_rsid_mask, "hg38_start"].astype(str)
                    + ":"
                    + variants_df.loc[missing_rsid_mask, "ref_allele"].astype(str)
                    + ":"
                    + variants_df.loc[missing_rsid_mask, "alt_allele"].astype(str)
                )
                snp_data.loc[missing_rsid_mask, "snp"] = coord_ids

            # Harmonize the SNPs
            harmonized_snp_data = harmonizer.harmonize_snp_batch(snp_data)

            if not harmonized_snp_data.empty:
                logger.info(
                    f"Successfully harmonized {len(harmonized_snp_data)} ClinVar SNPs"
                )
                return harmonized_snp_data
            else:
                logger.warning(
                    "Harmonization resulted in empty DataFrame, returning original data"
                )

        except Exception as e:
            logger.error(f"Error during ClinVar SNP harmonization: {e}")
            logger.info("Continuing with non-harmonized data")

    # If no harmonizer or harmonization failed, ensure consistent format
    # For variants without rsIDs, use coordinate-based format (NO clinvar_chr_pos)
    missing_rsid_mask = snp_data["snp"].isna()
    if missing_rsid_mask.any():
        coordinate_ids = (
            snp_data.loc[missing_rsid_mask, "hg38_chromosome"].astype(str)
            + ":"
            + snp_data.loc[missing_rsid_mask, "hg38_start"].astype(str)
            + ":"
            + variants_df.loc[missing_rsid_mask, "ref_allele"].astype(str)
            + ":"
            + variants_df.loc[missing_rsid_mask, "alt_allele"].astype(str)
        )
        snp_data.loc[missing_rsid_mask, "snp"] = coordinate_ids
        logger.info(
            f"Created coordinate-based IDs for {missing_rsid_mask.sum()} ClinVar variants without rsIDs"
        )

    return snp_data


def _extract_gene_symbol(geneinfo: str) -> Optional[str]:
    """
    Extract gene symbol from GENEINFO field.

    Args:
        geneinfo: GENEINFO field value (format: "GENE:ID|GENE:ID")

    Returns:
        First gene symbol or None
    """
    if not geneinfo or pd.isna(geneinfo):
        return None

    # GENEINFO format: "BRCA1:672|BRCA1:672"
    try:
        # Take first gene symbol before colon
        first_gene = geneinfo.split("|")[0].split(":")[0]
        return first_gene if first_gene else None
    except (IndexError, AttributeError):
        return None


def _standardize_chromosome(chrom: str) -> str:
    """
    Standardize chromosome notation.

    Args:
        chrom: Chromosome string

    Returns:
        Standardized chromosome (e.g., "chr1", "chrX")
    """
    if not chrom:
        return chrom

    # Remove 'chr' prefix if present, then add it back
    chrom_clean = chrom.replace("chr", "")

    # Handle special cases
    if chrom_clean in ["X", "Y", "M", "MT"]:
        return f"chr{chrom_clean}"

    # Handle numeric chromosomes
    try:
        chrom_num = int(chrom_clean)
        if 1 <= chrom_num <= 22:
            return f"chr{chrom_num}"
    except ValueError:
        pass

    # Return original if can't standardize
    return chrom
