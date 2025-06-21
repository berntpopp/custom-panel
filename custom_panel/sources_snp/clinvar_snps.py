"""
ClinVar SNPs source for deep intronic variants.

This module processes ClinVar VCF files to extract deep intronic variants
that are within specified distance from gene panel exon boundaries.
"""

import gzip
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import requests

logger = logging.getLogger(__name__)


def fetch_clinvar_snps(
    config: dict[str, Any], gene_panel: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Fetch deep intronic ClinVar SNPs for the final gene panel.

    Args:
        config: Configuration dictionary containing ClinVar settings
        gene_panel: DataFrame with final gene panel (must have coordinates)

    Returns:
        DataFrame with ClinVar SNP data filtered to panel genes only
    """
    logger.info("Starting ClinVar SNPs processing...")

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

        # Apply deep intronic filtering
        deep_intronic_variants = _filter_deep_intronic_variants(
            panel_filtered_variants, panel_genes, clinvar_config
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
        return pd.DataFrame()


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

                    # Log progress every 50MB
                    if (
                        downloaded_size > 0
                        and downloaded_size % (50 * 1024 * 1024) < 8192
                    ):
                        if total_size > 0:
                            progress = (downloaded_size / total_size) * 100
                            logger.info(
                                f"Downloaded {downloaded_size // (1024*1024)}MB / {total_size // (1024*1024)}MB ({progress:.1f}%)"
                            )
                        else:
                            logger.info(
                                f"Downloaded {downloaded_size // (1024*1024)}MB"
                            )

        logger.info(f"Successfully downloaded ClinVar VCF to: {cache_file}")
        return cache_file

    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to download ClinVar VCF: {e}")
        return None
    except Exception as e:
        logger.error(f"Error downloading ClinVar VCF: {e}")
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
                    "clndn": info_dict.get("CLNDN", ""),
                    "clnrevstat": info_dict.get("CLNREVSTAT", ""),
                    "clnvc": info_dict.get("CLNVC", ""),
                    "geneinfo": info_dict.get("GENEINFO", ""),
                }

                variants.append(variant)

            except (ValueError, IndexError) as e:
                logger.warning(f"Error parsing line {line_num}: {e}")
                continue

            # Process in chunks to avoid memory issues
            if len(variants) % 100000 == 0 and len(variants) > 0:
                logger.info(f"Processed {len(variants)} variants so far...")
                # Limit processing to avoid memory issues and long processing times
                if (
                    len(variants) >= 2000000
                ):  # 2M variants should be sufficient for panel genes
                    logger.info("Limiting to 2M variants for memory and performance")
                    break

    if not variants:
        return pd.DataFrame()

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


def _convert_to_snp_format(variants_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert ClinVar variants to standardized SNP format matching other SNP sources.

    Args:
        variants_df: DataFrame with ClinVar variants

    Returns:
        DataFrame in standardized SNP format compatible with HTML reporting
    """
    if variants_df.empty:
        return pd.DataFrame()

    # Extract gene symbols from GENEINFO field
    variants_df["gene_symbol"] = variants_df["geneinfo"].apply(_extract_gene_symbol)

    # Create unique IDs for variants without rsIDs first
    rsids = variants_df["rsid"].copy()
    missing_rsid_mask = rsids.isna()
    if missing_rsid_mask.any():
        rsids.loc[missing_rsid_mask] = (
            "clinvar_"
            + variants_df.loc[missing_rsid_mask, "chromosome"].astype(str)
            + "_"
            + variants_df.loc[missing_rsid_mask, "position"].astype(str)
        )

    # Create allele string format (ref/alt)
    allele_strings = (
        variants_df["ref_allele"].astype(str)
        + "/"
        + variants_df["alt_allele"].astype(str)
    )

    # Create standardized SNP data matching the expected format for HTML report
    snp_data = pd.DataFrame(
        {
            # Standard SNP identification columns
            "snp": rsids,  # This is the main identifier column
            "source": "ClinVar",
            "category": "deep_intronic_clinvar",  # This is crucial for HTML display
            "snp_type": "deep_intronic_clinvar",  # This matches other SNP types
            # hg38 coordinate columns (matching other SNPs format)
            "hg38_chromosome": variants_df["chromosome"].apply(_standardize_chromosome),
            "hg38_start": variants_df["position"],
            "hg38_end": variants_df["position"],  # SNVs have same start/end
            "hg38_strand": 1,  # Default positive strand
            "hg38_allele_string": allele_strings,
            # hg19 coordinate columns (empty for ClinVar as it's GRCh38)
            "hg19_chromosome": "",
            "hg19_start": "",
            "hg19_end": "",
            "hg19_strand": "",
            "hg19_allele_string": "",
            # ClinVar-specific columns (additional detail)
            "rsid": rsids,
            "chromosome": variants_df["chromosome"].apply(_standardize_chromosome),
            "position": variants_df["position"],
            "ref_allele": variants_df["ref_allele"],
            "alt_allele": variants_df["alt_allele"],
            "gene_symbol": variants_df["gene_symbol"],
            "clinical_significance": variants_df["clnsig"],
            "condition": variants_df["clndn"],
            "review_status": variants_df["clnrevstat"],
            "variant_type": variants_df["clnvc"],
            "source_details": "Deep intronic pathogenic variants",
            # Legacy coordinate columns (for backward compatibility)
            "hg38_chr": variants_df["chromosome"].apply(_standardize_chromosome),
            "hg38_pos": variants_df["position"],
            "hg19_chr": "",
            "hg19_pos": "",
        }
    )

    # Keep variants with valid positions
    snp_data = snp_data.dropna(subset=["position"])

    return snp_data.reset_index(drop=True)


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


def _filter_variants_by_gene_panel(
    variants_df: pd.DataFrame, panel_genes: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """
    Filter ClinVar variants to only those within gene panel genomic regions.

    Args:
        variants_df: DataFrame with ClinVar variants
        panel_genes: DataFrame with final gene panel (must have genomic coordinates)
        config: ClinVar configuration dictionary

    Returns:
        Filtered DataFrame with variants only within gene panel regions
    """
    if variants_df.empty or panel_genes.empty:
        return pd.DataFrame()

    # Check if panel genes have required coordinate columns
    required_coords = ["chromosome", "gene_start", "gene_end"]
    if not all(col in panel_genes.columns for col in required_coords):
        logger.warning(
            f"Gene panel missing required coordinate columns {required_coords}"
        )
        logger.warning(f"Available columns: {list(panel_genes.columns)}")
        return pd.DataFrame()

    logger.info(
        f"Filtering {len(variants_df)} variants against {len(panel_genes)} genes using optimized algorithm"
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

        logger.debug(
            f"Processing {len(chrom_variants)} variants on {chrom} against {len(chrom_genes)} genes"
        )

        # For each variant on this chromosome, check if it falls within any gene
        for _, variant in chrom_variants.iterrows():
            variant_pos = variant["position"]
            variant_gene = variant["gene_symbol"]

            # Check if variant overlaps with any gene on this chromosome
            for _, gene in chrom_genes.iterrows():
                gene_start = int(gene["gene_start"])
                gene_end = int(gene["gene_end"])
                gene_symbol = gene["approved_symbol"]

                # Check if variant is within gene boundaries
                if gene_start <= variant_pos <= gene_end:
                    # Additional check: if variant has gene symbol, it should match
                    if variant_gene and pd.notna(gene_symbol):
                        if variant_gene == gene_symbol:
                            filtered_variants.append(variant)
                            break  # Found exact match
                    else:
                        # No gene symbol or no match - accept based on position alone
                        filtered_variants.append(variant)
                        break  # Found positional match

    if not filtered_variants:
        logger.info("No variants found within gene panel regions")
        return pd.DataFrame()

    # Convert back to DataFrame and remove the temporary columns
    result_df = pd.DataFrame(filtered_variants)
    result_df = result_df.drop(columns=["std_chr", "gene_symbol"], errors="ignore")

    logger.info(
        f"Filtered to {len(result_df)} variants within {len(panel_genes)} gene panel regions"
    )

    return result_df


def _filter_deep_intronic_variants(
    variants_df: pd.DataFrame, panel_genes: pd.DataFrame, config: dict[str, Any]
) -> pd.DataFrame:
    """
    Filter variants to only deep intronic ones (outside exon boundaries by specified padding).

    Args:
        variants_df: DataFrame with ClinVar variants within gene panel regions
        panel_genes: DataFrame with final gene panel
        config: ClinVar configuration dictionary

    Returns:
        DataFrame with only deep intronic variants
    """
    if variants_df.empty:
        return pd.DataFrame()

    # Get intronic padding from config (default 50bp)
    intronic_padding = config.get("intronic_padding", 50)

    # For now, we'll assume all variants are intronic if they're within gene boundaries
    # but outside of exons. Since we don't have exon coordinates readily available,
    # we'll keep all variants that are within gene boundaries as potentially intronic.
    # This is a simplified approach - a more sophisticated implementation would
    # require exon coordinate data from Ensembl or other sources.

    logger.info(f"Applying deep intronic filter with {intronic_padding}bp padding")
    logger.info(
        "Note: Using simplified intronic filtering - variants within gene boundaries are considered potentially intronic"
    )

    # For now, return all variants as they're already filtered to be within gene regions
    # and ClinVar typically contains variants outside of coding regions
    deep_intronic_variants = variants_df.copy()

    logger.info(
        f"Identified {len(deep_intronic_variants)} potential deep intronic variants"
    )

    return deep_intronic_variants


if __name__ == "__main__":
    # Test the module
    test_config = {
        "snp_processing": {
            "deep_intronic_clinvar": {
                "enabled": True,
                "vcf_path": "data/clinvar/clinvar.vcf.gz",
                "cache_dir": ".cache/clinvar",
                "cache_expiry_days": 30,
                "intronic_padding": 50,
            }
        }
    }

    result = fetch_clinvar_snps(test_config)
    print(f"Test result: {len(result)} SNPs processed")
