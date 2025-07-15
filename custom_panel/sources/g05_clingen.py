"""
ClinGen data source extractor.

This module processes ClinGen gene validity evidence for cancer-related genes.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from ..core.cache_manager import CacheManager
from ..core.config_manager import ConfigManager
from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)


def _fetch_clingen_genes_api(
    url: str, cache_manager: CacheManager | None = None
) -> list[dict[str, Any]]:
    """
    Fetch ClinGen genes from the ClinGen API with caching support.

    Args:
        url: URL to the ClinGen API endpoint
        cache_manager: Optional cache manager for caching responses

    Returns:
        List of gene records with classification information

    Raises:
        ValueError: If the expected data structure is not found
        requests.RequestException: If the API request fails
    """
    # Check cache first if available
    if cache_manager and cache_manager.enabled:
        cached_data = cache_manager.get("clingen", url, "GET", None)
        if cached_data is not None:
            logger.info(f"Using cached ClinGen data from {url}")
            return cached_data

    logger.info(f"Attempting to fetch ClinGen genes from API: {url}")

    # Make request with User-Agent header
    headers = {
        "User-Agent": "custom-panel/1.0 (Python scientific tool for gene panel curation)",
        "Accept": "application/json",
    }

    try:
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
    except requests.RequestException as e:
        logger.error(f"Failed to fetch ClinGen API from {url}: {e}")
        raise

    try:
        # Parse JSON response
        data = response.json()

        # The API should return gene validity data
        gene_data = []

        # Handle different possible API response structures
        if isinstance(data, list):
            gene_data = data
        elif isinstance(data, dict):
            # Try common data keys
            for key in [
                "rows",
                "results",
                "data",
                "genes",
                "records",
                "affiliations",
                "geneValidities",
            ]:
                if key in data and isinstance(data[key], list):
                    gene_data = data[key]
                    break

        if not gene_data:
            raise ValueError("No gene data found in ClinGen API response")

        logger.info(f"Found {len(gene_data)} gene records from ClinGen API")

        # Filter out refuted classifications
        filtered_genes = []
        for record in gene_data:
            if not isinstance(record, dict):
                continue

            classification = record.get("classification", "").strip()
            if classification and classification.lower() != "refuted":
                filtered_genes.append(record)

        logger.info(
            f"Filtered out refuted genes: {len(filtered_genes)} remaining from {len(gene_data)} total"
        )

        # Cache the processed data
        if cache_manager and cache_manager.enabled:
            cache_manager.set("clingen", url, "GET", None, filtered_genes)
            logger.info("Cached ClinGen data for future use")

        return filtered_genes

    except (json.JSONDecodeError, KeyError, TypeError) as e:
        raise ValueError(f"Failed to parse ClinGen API response: {e}") from e


def fetch_clingen_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch ClinGen gene validity data with caching support.

    Prioritizes live scraping from ClinGen website, with fallback to local file.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with ClinGen gene validity data
    """
    config_manager = ConfigManager(config)
    clingen_config = config_manager.get_source_config("ClinGen")

    if not clingen_config.get("enabled", True):
        logger.info("ClinGen data source is disabled")
        return pd.DataFrame()

    # Initialize cache manager
    cache_config = config.get("performance", {})
    cache_manager = CacheManager(
        cache_dir=cache_config.get("cache_dir", ".cache"),
        cache_ttl=cache_config.get("cache_ttl", 2592000),  # 30 days
        enabled=cache_config.get("enable_caching", True),
    )

    base_evidence_score = clingen_config.get("evidence_score", 1.2)
    classification_scores = clingen_config.get(
        "classification_scores",
        {
            "Definitive": 2.0,
            "Strong": 1.5,
            "Moderate": 1.0,
            "Limited": 0.5,
            "Disputed": 0.3,
            "Refuted": 0.0,
        },
    )
    url = clingen_config.get("url")
    file_path = clingen_config.get("file_path")
    retrieval_date = datetime.now().strftime("%Y-%m-%d")

    genes = []
    evidence_scores = []
    source_details = []

    # Priority 1: Try to fetch from live API with caching
    if url:
        try:
            gene_records = _fetch_clingen_genes_api(url, cache_manager)

            # Extract gene symbols and classifications
            for record in gene_records:
                # Look for gene symbol in various possible fields
                gene_symbol = None
                for field in [
                    "symbol",
                    "gene_symbol",
                    "gene",
                    "hgnc_symbol",
                    "approved_symbol",
                    "geneSymbol",
                    "hgncSymbol",
                ]:
                    if field in record and record[field]:
                        gene_symbol = str(record[field]).strip()
                        break

                if not gene_symbol:
                    continue

                # Look for classification information
                classification = "Unknown"
                for field in [
                    "classification",
                    "classification_title",
                    "validity",
                    "evidence_level",
                ]:
                    if field in record and record[field]:
                        classification = str(record[field]).strip()
                        break

                # Calculate evidence score based on classification
                classification_multiplier = classification_scores.get(
                    classification, 1.0
                )
                evidence_score = base_evidence_score * classification_multiplier

                genes.append(gene_symbol)
                evidence_scores.append(evidence_score)
                source_details.append(
                    f"URL:{url}|Date:{retrieval_date}|Classification:{classification}"
                )

            if genes:
                logger.info(
                    f"Successfully fetched {len(genes)} ClinGen genes from live URL: {url}"
                )

        except Exception as e:
            logger.warning(f"Failed to fetch ClinGen genes from API {url}: {e}")
            logger.warning("Attempting to use fallback data source")

    # Priority 2: Fall back to local file if scraping failed
    if not genes and file_path:
        try:
            if Path(file_path).exists():
                file_path_obj = Path(file_path)

                if file_path_obj.suffix.lower() in [".csv", ".tsv"]:
                    # Handle CSV/TSV files
                    separator = "\t" if file_path_obj.suffix.lower() == ".tsv" else ","
                    df_fallback = pd.read_csv(file_path, sep=separator)

                    # Try common column names for gene symbols
                    gene_columns = [
                        "gene_symbol",
                        "Gene Symbol",
                        "Gene",
                        "symbol",
                        "approved_symbol",
                        "hgnc_symbol",
                        "HGNC Symbol",
                    ]
                    gene_column = None
                    for col in gene_columns:
                        if col in df_fallback.columns:
                            gene_column = col
                            break

                    if gene_column:
                        # Try to find classification column
                        classification_columns = [
                            "classification",
                            "Classification",
                            "validity",
                            "Validity",
                            "evidence_level",
                            "Evidence Level",
                            "category",
                            "Category",
                        ]
                        classification_column = None
                        for col in classification_columns:
                            if col in df_fallback.columns:
                                classification_column = col
                                break

                        for _, row in df_fallback.iterrows():
                            gene_symbol = str(row[gene_column]).strip()
                            if gene_symbol and gene_symbol not in ["", "nan", "NaN"]:
                                classification = "Unknown"
                                if classification_column and pd.notna(
                                    row[classification_column]
                                ):
                                    classification = str(
                                        row[classification_column]
                                    ).strip()

                                # Skip refuted genes
                                if classification.lower() == "refuted":
                                    continue

                                # Calculate evidence score based on classification
                                classification_multiplier = classification_scores.get(
                                    classification, 1.0
                                )
                                evidence_score = (
                                    base_evidence_score * classification_multiplier
                                )

                                genes.append(gene_symbol)
                                evidence_scores.append(evidence_score)
                                source_details.append(
                                    f"File:{file_path_obj.name}|Classification:{classification}"
                                )

                        logger.info(
                            f"Loaded {len(genes)} ClinGen genes from fallback file: {file_path}"
                        )
                    else:
                        logger.warning(f"No suitable gene column found in {file_path}")

                elif file_path_obj.suffix.lower() in [".xlsx", ".xls"]:
                    # Handle Excel files
                    df_fallback = pd.read_excel(file_path)

                    # Try common column names for gene symbols
                    gene_columns = [
                        "gene_symbol",
                        "Gene Symbol",
                        "Gene",
                        "symbol",
                        "approved_symbol",
                        "hgnc_symbol",
                        "HGNC Symbol",
                    ]
                    gene_column = None
                    for col in gene_columns:
                        if col in df_fallback.columns:
                            gene_column = col
                            break

                    if gene_column:
                        # Try to find classification column
                        classification_columns = [
                            "classification",
                            "Classification",
                            "validity",
                            "Validity",
                            "evidence_level",
                            "Evidence Level",
                            "category",
                            "Category",
                        ]
                        classification_column = None
                        for col in classification_columns:
                            if col in df_fallback.columns:
                                classification_column = col
                                break

                        for _, row in df_fallback.iterrows():
                            gene_symbol = str(row[gene_column]).strip()
                            if gene_symbol and gene_symbol not in ["", "nan", "NaN"]:
                                classification = "Unknown"
                                if classification_column and pd.notna(
                                    row[classification_column]
                                ):
                                    classification = str(
                                        row[classification_column]
                                    ).strip()

                                # Skip refuted genes
                                if classification.lower() == "refuted":
                                    continue

                                # Calculate evidence score based on classification
                                classification_multiplier = classification_scores.get(
                                    classification, 1.0
                                )
                                evidence_score = (
                                    base_evidence_score * classification_multiplier
                                )

                                genes.append(gene_symbol)
                                evidence_scores.append(evidence_score)
                                source_details.append(
                                    f"File:{file_path_obj.name}|Classification:{classification}"
                                )

                        logger.info(
                            f"Loaded {len(genes)} ClinGen genes from fallback file: {file_path}"
                        )
                    else:
                        logger.warning(f"No suitable gene column found in {file_path}")

        except Exception as e:
            logger.warning(
                f"Could not load ClinGen genes from fallback file {file_path}: {e}"
            )

    if not genes:
        logger.error("No ClinGen genes available from any source")
        return pd.DataFrame()

    # Remove duplicates while preserving highest evidence score for each gene
    gene_dict: dict[str, tuple[float, str]] = {}

    for gene, score, detail in zip(
        genes, evidence_scores, source_details, strict=False
    ):
        if gene not in gene_dict or score > gene_dict[gene][0]:
            gene_dict[gene] = (score, detail)

    # Extract unique genes with their best scores
    unique_genes = list(gene_dict.keys())
    unique_scores = [gene_dict[gene][0] for gene in unique_genes]
    unique_details = [gene_dict[gene][1] for gene in unique_genes]

    # Create standardized dataframe
    df = create_standard_dataframe(
        genes=unique_genes,
        source_name="ClinGen",
        evidence_scores=unique_scores,
        source_details=unique_details,
        gene_names_reported=unique_genes,
    )

    logger.info(f"Created ClinGen dataset with {len(df)} genes")
    return df
