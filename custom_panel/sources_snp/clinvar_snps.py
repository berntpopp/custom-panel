"""
Optimized ClinVar SNPs source using tabix-based streaming and existing exon data.

This module processes ClinVar VCF files to extract deep intronic variants
using efficient tabix queries and reusing exon data from the annotation phase.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import requests

try:
    import pysam
except ImportError:
    pysam = None  # type: ignore

logger = logging.getLogger(__name__)


class IntervalTree:
    """Simple interval tree implementation for efficient overlap detection."""

    def __init__(self) -> None:
        self.intervals: dict[str, list[tuple[int, int]]] = {}

    def add(self, chromosome: str, start: int, end: int) -> None:
        """Add an interval to the tree."""
        if chromosome not in self.intervals:
            self.intervals[chromosome] = []
        self.intervals[chromosome].append((start, end))

    def sort(self) -> None:
        """Sort intervals for efficient querying."""
        for chromosome in self.intervals:
            self.intervals[chromosome].sort()

    def overlaps(self, chromosome: str, position: int) -> bool:
        """Check if position overlaps with any interval on chromosome."""
        if chromosome not in self.intervals:
            return False

        # Binary search would be more efficient, but linear is sufficient for our use case
        for start, end in self.intervals[chromosome]:
            if start <= position <= end:
                return True
        return False


def fetch_clinvar_snps(
    config: dict[str, Any],
    gene_panel: Optional[pd.DataFrame] = None,
    ensembl_client: Optional[Any] = None,
    harmonizer: Optional[Any] = None,
    transcript_data: Optional[dict[str, Any]] = None,
) -> pd.DataFrame:
    """
    Fetch deep intronic ClinVar SNPs using optimized tabix-based streaming.

    Args:
        config: Configuration dictionary containing ClinVar settings
        gene_panel: DataFrame with final gene panel (must have coordinates)
        ensembl_client: Optional Ensembl client (not used in optimized version)
        harmonizer: Optional SNP harmonizer for coordinate conversion
        transcript_data: Pre-computed transcript data from annotation phase

    Returns:
        DataFrame with ClinVar SNP data filtered to deep intronic variants only
    """
    logger.info("Starting optimized ClinVar SNPs processing using tabix streaming...")

    if pysam is None:
        logger.error(
            "pysam is required for tabix support. Install with: pip install pysam"
        )
        return pd.DataFrame()

    clinvar_config = config.get("snp_processing", {}).get("deep_intronic_clinvar", {})

    if not clinvar_config.get("enabled", False):
        logger.info("ClinVar SNPs processing is disabled")
        return pd.DataFrame()

    # Check if we have gene panel data and transcript data
    if gene_panel is None or gene_panel.empty:
        logger.warning("No gene panel provided for ClinVar filtering")
        return pd.DataFrame()

    if transcript_data is None:
        logger.warning("No transcript data provided for exon region extraction")
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
        return pd.DataFrame()

    # Remove genes without coordinate data
    panel_genes = panel_genes.dropna(subset=required_coords)
    if panel_genes.empty:
        logger.warning("No genes with valid coordinate data")
        return pd.DataFrame()

    try:
        # Step 1: Extract exon regions from existing transcript data
        intronic_padding = clinvar_config.get("intronic_padding", 50)
        exon_intervals = _extract_exon_regions_from_transcript_data(
            transcript_data, panel_genes, intronic_padding
        )

        if not exon_intervals.intervals:
            logger.warning("No exon regions extracted from transcript data")
            return pd.DataFrame()

        # Step 2: Get or download ClinVar VCF file with tabix index
        vcf_file, tbi_file = _get_clinvar_vcf_with_index(clinvar_config)
        if not vcf_file or not vcf_file.exists():
            logger.error("Could not obtain ClinVar VCF file")
            return pd.DataFrame()

        if not tbi_file or not tbi_file.exists():
            logger.error("Could not obtain ClinVar tabix index file")
            return pd.DataFrame()

        # Step 3: Stream ClinVar variants by gene regions using tabix
        logger.info("Streaming ClinVar variants using tabix queries...")
        variants = list(
            _stream_clinvar_variants_by_genes(vcf_file, panel_genes, clinvar_config)
        )

        if not variants:
            logger.warning("No variants found in ClinVar for gene panel regions")
            return pd.DataFrame()

        logger.info(f"Found {len(variants)} variants in gene panel regions")

        # Step 4: Filter for deep intronic variants (outside padded exon regions)
        logger.info("Filtering for deep intronic variants...")
        deep_intronic_variants = _filter_deep_intronic_streaming(
            variants, exon_intervals
        )

        if not deep_intronic_variants:
            logger.warning("No deep intronic variants found after filtering")
            return pd.DataFrame()

        logger.info(
            f"Found {len(deep_intronic_variants)} deep intronic ClinVar variants"
        )

        # Step 5: Convert to standardized SNP format
        snp_formatted = _convert_to_snp_format(deep_intronic_variants, harmonizer)

        logger.info(f"Successfully processed {len(snp_formatted)} ClinVar SNPs")
        return snp_formatted

    except Exception as e:
        logger.error(f"Error processing ClinVar SNPs: {str(e)}")
        return pd.DataFrame()


def _extract_exon_regions_from_transcript_data(
    transcript_data: dict[str, Any],
    panel_genes: pd.DataFrame,
    intronic_padding: int,
) -> IntervalTree:
    """
    Extract padded exon regions from existing transcript data.

    Args:
        transcript_data: Pre-computed transcript data from annotation phase
        panel_genes: DataFrame with panel genes
        intronic_padding: Padding around exons in base pairs

    Returns:
        IntervalTree with padded exon regions
    """
    logger.info(
        f"Extracting exon regions with {intronic_padding}bp padding from transcript data"
    )

    exon_intervals = IntervalTree()
    genes_with_exons = 0
    total_exons = 0

    for _, gene_row in panel_genes.iterrows():
        gene_symbol = gene_row["approved_symbol"]

        # Get transcript data for this gene
        gene_data = transcript_data.get(gene_symbol)
        if not gene_data or "all_transcripts" not in gene_data:
            continue

        # Process canonical transcript if available
        canonical_transcript_id = gene_row.get("canonical_transcript")
        if pd.notna(canonical_transcript_id):
            transcript = _find_transcript_by_id(
                gene_data["all_transcripts"], canonical_transcript_id
            )
            if transcript and "Exon" in transcript:
                exons = _extract_exons_from_transcript(
                    transcript, gene_symbol, intronic_padding
                )

                for exon in exons:
                    exon_intervals.add(exon["chromosome"], exon["start"], exon["end"])
                    total_exons += 1

                if exons:
                    genes_with_exons += 1

    # Sort intervals for efficient querying
    exon_intervals.sort()

    logger.info(f"Extracted {total_exons} exon regions from {genes_with_exons} genes")
    return exon_intervals


def _find_transcript_by_id(
    transcripts: list[dict[str, Any]], transcript_id: str
) -> Optional[dict[str, Any]]:
    """Find transcript by ID in the transcripts list."""
    for transcript in transcripts:
        if transcript.get("id") == transcript_id:
            return transcript
    return None


def _extract_exons_from_transcript(
    transcript: dict[str, Any], gene_symbol: str, intronic_padding: int
) -> list[dict[str, Any]]:
    """Extract exon coordinates from transcript data."""
    exons: list[dict[str, Any]] = []

    if "Exon" not in transcript:
        return exons

    for exon in transcript["Exon"]:
        if all(key in exon for key in ["seq_region_name", "start", "end"]):
            # Apply padding to create exclusion zones
            padded_start = max(1, exon["start"] - intronic_padding)
            padded_end = exon["end"] + intronic_padding

            exon_info = {
                "chromosome": _standardize_chromosome(exon["seq_region_name"]),
                "start": padded_start,
                "end": padded_end,
                "gene_symbol": gene_symbol,
                "exon_id": exon.get("id", ""),
                "rank": exon.get("rank", 0),
            }
            exons.append(exon_info)

    return exons


def _validate_vcf_file(vcf_file: Path, tbi_file: Path) -> bool:
    """
    Validate that VCF file and index are not corrupted.

    Args:
        vcf_file: Path to VCF file
        tbi_file: Path to tabix index file

    Returns:
        True if files are valid, False otherwise
    """
    if not vcf_file.exists() or not tbi_file.exists():
        return False

    try:
        # Try to open the tabix file and read a few records
        tbx = pysam.TabixFile(str(vcf_file))

        # Try to iterate through chromosomes to test file integrity
        for chrom in ["1", "2", "3"]:  # Test a few chromosomes
            try:
                # Just try to create an iterator - don't need to read all records
                records = tbx.fetch(chrom, 1, 1000000)
                # Try to read first record to test if file is readable
                next(records, None)
                break  # If we can read from any chromosome, file is likely okay
            except StopIteration:
                # No records in this region, try next chromosome
                continue
            except Exception:
                # If we can't read from any chromosome, file is corrupted
                return False

        tbx.close()
        return True

    except Exception as e:
        logger.warning(f"VCF file validation failed: {e}")
        return False


def _get_clinvar_vcf_with_index(
    config: dict[str, Any],
) -> tuple[Optional[Path], Optional[Path]]:
    """
    Get ClinVar VCF file and tabix index - either from local path or by downloading.

    Args:
        config: ClinVar configuration dictionary

    Returns:
        Tuple of (VCF file path, tabix index file path) or (None, None) if failed
    """
    # Check if local file path is provided
    vcf_path = config.get("vcf_path")
    if vcf_path:
        vcf_file = Path(vcf_path)
        tbi_file = Path(str(vcf_path) + ".tbi")

        if _validate_vcf_file(vcf_file, tbi_file):
            logger.info(f"Using provided ClinVar VCF file: {vcf_path}")
            return vcf_file, tbi_file
        else:
            logger.warning(
                f"Provided ClinVar VCF or index file not found or corrupted: {vcf_path}"
            )

    # Download ClinVar VCF and index if not provided or not found
    vcf_url = config.get("vcf_url")
    if not vcf_url:
        logger.error("No ClinVar VCF URL provided for download")
        return None, None

    cache_dir = Path(config.get("cache_dir", ".cache/clinvar"))
    cache_expiry_days = config.get("cache_expiry_days", 30)

    # Create cache directory
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Download VCF file and index
    downloaded_vcf = _download_clinvar_vcf(vcf_url, cache_dir, cache_expiry_days)
    downloaded_tbi = _download_clinvar_index(vcf_url, cache_dir, cache_expiry_days)

    # Validate downloaded files
    if downloaded_vcf and downloaded_tbi:
        if _validate_vcf_file(downloaded_vcf, downloaded_tbi):
            return downloaded_vcf, downloaded_tbi
        else:
            logger.warning(
                "Downloaded ClinVar files are corrupted, removing cache and retrying..."
            )
            # Remove corrupted files
            if downloaded_vcf.exists():
                downloaded_vcf.unlink()
            if downloaded_tbi.exists():
                downloaded_tbi.unlink()

            # Retry download once
            downloaded_vcf = _download_clinvar_vcf(
                vcf_url, cache_dir, cache_expiry_days
            )
            downloaded_tbi = _download_clinvar_index(
                vcf_url, cache_dir, cache_expiry_days
            )

            if (
                downloaded_vcf
                and downloaded_tbi
                and _validate_vcf_file(downloaded_vcf, downloaded_tbi)
            ):
                return downloaded_vcf, downloaded_tbi

    return None, None


def _download_clinvar_vcf(
    url: str, cache_dir: Path, cache_expiry_days: int
) -> Optional[Path]:
    """Download ClinVar VCF file with caching."""
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

    # Download the file
    logger.info(f"Downloading ClinVar VCF from: {url}")
    try:
        response = requests.get(url, stream=True, timeout=300)
        response.raise_for_status()

        with open(cache_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        logger.info(f"Successfully downloaded ClinVar VCF to: {cache_file}")
        return cache_file

    except Exception as e:
        logger.error(f"Error downloading ClinVar VCF: {str(e)}")
        return None


def _download_clinvar_index(
    vcf_url: str, cache_dir: Path, cache_expiry_days: int
) -> Optional[Path]:
    """Download ClinVar tabix index file with caching."""
    # Construct index URL
    index_url = vcf_url + ".tbi"
    index_filename = vcf_url.split("/")[-1] + ".tbi"

    cache_file = cache_dir / index_filename

    # Check if cached file exists and is recent enough
    if cache_file.exists():
        file_age = datetime.now() - datetime.fromtimestamp(cache_file.stat().st_mtime)
        if file_age.days < cache_expiry_days:
            logger.info(f"Using cached ClinVar index: {cache_file}")
            return cache_file

    # Download the index file
    logger.info(f"Downloading ClinVar tabix index from: {index_url}")
    try:
        response = requests.get(index_url, stream=True, timeout=300)
        response.raise_for_status()

        with open(cache_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        logger.info(f"Successfully downloaded ClinVar index to: {cache_file}")
        return cache_file

    except Exception as e:
        logger.error(f"Error downloading ClinVar index: {str(e)}")
        return None


def _stream_clinvar_variants_by_genes(
    vcf_file: Path, panel_genes: pd.DataFrame, config: dict[str, Any]
) -> Iterator[dict[str, Any]]:
    """
    Stream ClinVar variants for gene panel regions using tabix.

    Args:
        vcf_file: Path to ClinVar VCF file
        panel_genes: DataFrame with gene panel regions
        config: ClinVar configuration

    Yields:
        Variant dictionaries for pathogenic variants in gene regions
    """
    max_indel_size = config.get("max_indel_size", 50)

    try:
        # Open tabix file
        tbx = pysam.TabixFile(str(vcf_file))

        variants_processed = 0
        variants_yielded = 0

        # Query each gene region
        for _, gene in panel_genes.iterrows():
            chromosome = _standardize_chromosome_for_clinvar(str(gene["chromosome"]))
            start = int(gene["gene_start"])
            end = int(gene["gene_end"])

            try:
                # Query tabix for this genomic region (convert to 0-based coordinates)
                for record in tbx.fetch(chromosome, start - 1, end):
                    variants_processed += 1

                    # Parse VCF record
                    variant = _parse_vcf_record(record)
                    if not variant:
                        continue

                    # Apply pathogenic filtering during streaming
                    if not _is_pathogenic_variant(variant):
                        continue

                    # Apply size filtering during streaming
                    if not _is_small_variant(variant, max_indel_size):
                        continue

                    variants_yielded += 1
                    yield variant

            except Exception as e:
                logger.warning(f"Error querying region {chromosome}:{start}-{end}: {e}")
                continue

        logger.info(
            f"Processed {variants_processed} variants, yielded {variants_yielded} pathogenic variants"
        )

    except Exception as e:
        logger.error(f"Error opening tabix file {vcf_file}: {e}")
        return


def _parse_vcf_record(record: str) -> Optional[dict[str, Any]]:
    """Parse a VCF record string into a variant dictionary."""
    try:
        fields = record.strip().split("\t")
        if len(fields) < 8:
            return None

        chrom = fields[0]
        pos = int(fields[1])
        variant_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        info = fields[7]

        # Parse INFO field
        info_dict = _parse_info_field(info)

        return {
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

    except (ValueError, IndexError) as e:
        logger.warning(f"Error parsing VCF record: {str(e)}")
        return None


def _parse_info_field(info: str) -> dict[str, str]:
    """Parse VCF INFO field into dictionary."""
    info_dict = {}

    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = "True"

    return info_dict


def _is_pathogenic_variant(variant: dict[str, Any]) -> bool:
    """Check if variant is pathogenic or likely pathogenic."""
    clnsig = variant.get("clnsig", "")

    pathogenic_values = [
        "Pathogenic",
        "Likely_pathogenic",
        "Pathogenic/Likely_pathogenic",
        "Likely_pathogenic/Pathogenic",
    ]

    # Check for exact matches
    if clnsig in pathogenic_values:
        return True

    # Check for pathogenic keywords in complex classifications
    if "Pathogenic" in clnsig or "Likely_pathogenic" in clnsig:
        return True

    return False


def _is_small_variant(variant: dict[str, Any], max_indel_size: int) -> bool:
    """Check if variant is small enough (not a large CNV)."""
    ref_len = len(variant.get("ref_allele", ""))
    alt_len = len(variant.get("alt_allele", ""))
    max_len = max(ref_len, alt_len)

    return max_len <= max_indel_size


def _filter_deep_intronic_streaming(
    variants: list[dict[str, Any]], exon_intervals: IntervalTree
) -> list[dict[str, Any]]:
    """
    Filter variants to only deep intronic ones using interval tree.

    Args:
        variants: List of variant dictionaries
        exon_intervals: IntervalTree with padded exon regions

    Returns:
        List of deep intronic variants (outside padded exon regions)
    """
    deep_intronic_variants = []

    for variant in variants:
        chromosome = _standardize_chromosome(variant["chromosome"])
        position = variant["position"]

        # Check if variant overlaps with any padded exon region
        if not exon_intervals.overlaps(chromosome, position):
            # Variant is outside all padded exon regions = deep intronic
            deep_intronic_variants.append(variant)

    return deep_intronic_variants


def _convert_to_snp_format(
    variants: list[dict[str, Any]], harmonizer: Optional[Any] = None
) -> pd.DataFrame:
    """Convert ClinVar variants to standardized SNP format."""
    if not variants:
        return pd.DataFrame()

    # Convert to DataFrame for easier processing
    variants_df = pd.DataFrame(variants)

    # Extract gene symbols from GENEINFO field
    variants_df["gene_symbol"] = variants_df["geneinfo"].apply(_extract_gene_symbol)

    # Create allele string format
    allele_strings = (
        variants_df["ref_allele"].astype(str)
        + "/"
        + variants_df["alt_allele"].astype(str)
    )

    # Create initial SNP data with hg38 coordinates
    snp_data = pd.DataFrame(
        {
            "snp": variants_df["rsid"],
            "source": "ClinVar",
            "category": "deep_intronic_clinvar",
            "snp_type": "deep_intronic_clinvar",
            "hg38_chromosome": variants_df["chromosome"].apply(_standardize_chromosome),
            "hg38_start": variants_df["position"],
            "hg38_end": variants_df["position"],
            "hg38_strand": 1,
            "hg38_allele_string": allele_strings,
            # Only hg38 coordinates are supported
            "clinical_significance": variants_df["clnsig"],
            "gene_symbol": variants_df["gene_symbol"],
            "condition": variants_df["clndn"],
            "review_status": variants_df["clnrevstat"],
            "variant_type": variants_df["clnvc"],
            "origin": variants_df["origin"],
            "molecular_consequence": variants_df["clnvcso"],
            "variant_id": variants_df["clnvi"],
            "processing_date": datetime.now().strftime("%Y-%m-%d"),
            "source_file": "ClinVar VCF",
        }
    )

    # Apply harmonization if available
    if harmonizer is not None:
        try:
            logger.info(f"Harmonizing {len(snp_data)} ClinVar SNPs")

            # Handle missing rsIDs
            missing_rsid_mask = snp_data["snp"].isna()
            if missing_rsid_mask.any():
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

        except Exception as e:
            logger.error(f"Error during ClinVar SNP harmonization: {e}")

    # Handle missing rsIDs for non-harmonized data
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

    return snp_data


def _extract_gene_symbol(geneinfo: str) -> Optional[str]:
    """Extract gene symbol from GENEINFO field."""
    if not geneinfo or pd.isna(geneinfo):
        return None

    try:
        # GENEINFO format: "BRCA1:672|BRCA1:672"
        first_gene = geneinfo.split("|")[0].split(":")[0]
        return first_gene if first_gene else None
    except (IndexError, AttributeError):
        return None


def _standardize_chromosome(chrom: str) -> str:
    """Standardize chromosome notation to chr prefix format."""
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


def _standardize_chromosome_for_clinvar(chrom: str) -> str:
    """Standardize chromosome notation for ClinVar VCF format (no chr prefix)."""
    if not chrom:
        return chrom

    # Remove 'chr' prefix if present
    chrom_clean = chrom.replace("chr", "")

    # Handle special cases
    if chrom_clean in ["X", "Y", "M", "MT"]:
        return chrom_clean

    # Handle numeric chromosomes
    try:
        chrom_num = int(chrom_clean)
        if 1 <= chrom_num <= 22:
            return str(chrom_num)
    except ValueError:
        pass

    # Return cleaned version if can't standardize
    return chrom_clean
