"""
GenCC (Gene Curation Coalition) data source extractor.

This module processes GenCC gene-disease associations for cancer-related genes.
"""

import logging
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)


def fetch_gencc_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch GenCC gene-disease association data.

    Prioritizes live download from GenCC website, with fallback to local file.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with GenCC gene-disease data
    """
    gencc_config = config.get("data_sources", {}).get("TheGenCC", {})

    if not gencc_config.get("enabled", True):
        logger.info("TheGenCC data source is disabled")
        return pd.DataFrame()

    base_evidence_score = gencc_config.get("evidence_score", 1.2)
    classification_scores = gencc_config.get(
        "classification_scores",
        {
            "Definitive": 2.0,
            "Strong": 1.5,
            "Moderate": 1.0,
            "Supportive": 0.7,
            "Limited": 0.5,
            "Disputed Evidence": 0.3,
            "No Known Disease Relationship": 0.0,
            "Refuted Evidence": 0.0,
        },
    )
    url = gencc_config.get("url")
    file_path = gencc_config.get("file_path")
    filter_keywords = gencc_config.get("filter_keywords", ["cancer", "tumor", "tumour"])
    retrieval_date = datetime.now().strftime("%Y-%m-%d")

    df_source = None
    source_prefix = ""

    # Priority 1: Try to download from live URL
    if url:
        try:
            logger.info(f"Attempting to download GenCC data from: {url}")

            # Make request with User-Agent header
            headers = {
                "User-Agent": "custom-panel/1.0 (Python scientific tool for gene panel curation)"
            }

            response = requests.get(url, headers=headers, timeout=60)
            response.raise_for_status()

            # Read Excel data directly from response content
            df_source = pd.read_excel(BytesIO(response.content), engine="openpyxl")
            source_prefix = f"URL:{url}|Date:{retrieval_date}"

            logger.info(
                f"Successfully downloaded GenCC data with {len(df_source)} records from live URL"
            )

        except Exception as e:
            logger.warning(f"Failed to download GenCC data from {url}: {e}")
            logger.warning("Attempting to use fallback data source")

    # Priority 2: Fall back to local file if download failed
    if df_source is None and file_path:
        try:
            if Path(file_path).exists():
                file_path_obj = Path(file_path)

                if file_path_obj.suffix.lower() in [".xlsx", ".xls"]:
                    df_source = pd.read_excel(file_path, engine="openpyxl")
                elif file_path_obj.suffix.lower() == ".csv":
                    df_source = pd.read_csv(file_path)
                elif file_path_obj.suffix.lower() == ".tsv":
                    df_source = pd.read_csv(file_path, sep="\t")
                else:
                    raise ValueError(f"Unsupported file format: {file_path_obj.suffix}")

                source_prefix = f"File:{file_path_obj.name}"
                logger.info(
                    f"Loaded GenCC data with {len(df_source)} records from fallback file: {file_path}"
                )

        except Exception as e:
            logger.warning(
                f"Could not load GenCC data from fallback file {file_path}: {e}"
            )

    if df_source is None:
        logger.error("No GenCC data available from any source")
        return pd.DataFrame()

    # Filter for cancer-related diseases
    logger.info(
        f"Filtering GenCC data for cancer-related diseases using keywords: {filter_keywords}"
    )

    # Look for disease information in various possible columns
    disease_columns = [
        "disease_original_title",
        "disease_title",
        "disease",
        "condition",
        "phenotype",
        "disease_name",
        "disease_label",
    ]

    cancer_mask = pd.Series([False] * len(df_source))

    for col in disease_columns:
        if col in df_source.columns:
            # Case-insensitive search for cancer keywords
            for keyword in filter_keywords:
                mask = (
                    df_source[col]
                    .astype(str)
                    .str.contains(keyword, case=False, na=False)
                )
                cancer_mask = cancer_mask | mask
            break  # Use the first available disease column

    if not cancer_mask.any():
        logger.warning("No cancer-related diseases found in GenCC data")
        return pd.DataFrame()

    df_filtered = df_source[cancer_mask].copy()
    logger.info(f"Filtered to {len(df_filtered)} cancer-related records")

    # Extract gene symbols
    gene_columns = [
        "gene_symbol",
        "Gene Symbol",
        "gene",
        "symbol",
        "hgnc_symbol",
        "HGNC Symbol",
        "approved_symbol",
    ]

    gene_column = None
    for col in gene_columns:
        if col in df_filtered.columns:
            gene_column = col
            break

    if not gene_column:
        logger.error(
            f"No suitable gene column found in GenCC data. Available columns: {list(df_filtered.columns)}"
        )
        return pd.DataFrame()

    # Try to find classification column
    classification_columns = [
        "classification_title",
        "classification",
        "evidence_level",
        "validity",
        "category",
        "confidence",
    ]

    classification_column = None
    for col in classification_columns:
        if col in df_filtered.columns:
            classification_column = col
            break

    # Extract unique genes and their details
    genes = []
    evidence_scores = []
    source_details = []

    for _, row in df_filtered.iterrows():
        gene_symbol = str(row[gene_column]).strip()

        # Skip empty or invalid gene symbols
        if not gene_symbol or gene_symbol in ["", "nan", "NaN", "None"]:
            continue

        # Get classification if available
        classification = "Unknown"
        if classification_column and pd.notna(row[classification_column]):
            classification = str(row[classification_column]).strip()

        # Filter out low-quality classifications
        classification_multiplier = classification_scores.get(classification, 1.0)
        if classification_multiplier == 0.0:
            continue  # Skip "No Known Disease Relationship" and "Refuted Evidence"

        # Calculate evidence score based on classification
        evidence_score = base_evidence_score * classification_multiplier

        genes.append(gene_symbol)
        evidence_scores.append(evidence_score)
        source_details.append(f"{source_prefix}|Classification:{classification}")

    if not genes:
        logger.error("No valid gene symbols found in filtered GenCC data")
        return pd.DataFrame()

    # Remove duplicates while preserving highest evidence score for each gene
    gene_dict = {}

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
        source_name="TheGenCC",
        evidence_scores=unique_scores,
        source_details=unique_details,
        gene_names_reported=unique_genes,
    )

    logger.info(f"Created GenCC dataset with {len(df)} unique genes")
    return df
