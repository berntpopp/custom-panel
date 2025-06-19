"""
ACMG Incidental Findings data source extractor.

This module processes ACMG recommendations for reporting incidental findings
in clinical exome and genome sequencing.
"""

import logging
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from bs4 import BeautifulSoup, Tag

from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)

# ACMG SF v3.2 genes (as of 2023) - can be overridden by config file
DEFAULT_ACMG_GENES = [
    "BRCA1",
    "BRCA2",
    "TP53",
    "STK11",
    "MLH1",
    "MSH2",
    "MSH6",
    "PMS2",
    "EPCAM",
    "VHL",
    "MEN1",
    "RET",
    "PTEN",
    "TSC1",
    "TSC2",
    "WT1",
    "NF2",
    "COL3A1",
    "FBN1",
    "TGFBR1",
    "TGFBR2",
    "SMAD3",
    "ACTA2",
    "MYLK",
    "MYH11",
    "KCNQ1",
    "KCNH2",
    "SCN5A",
    "LDLR",
    "APOB",
    "PCSK9",
    "RYR1",
    "CACNA1S",
    "ATP7B",
    "BMPR1A",
    "SMAD4",
    "MUTYH",
    "SDHD",
    "SDHAF2",
    "SDHC",
    "SDHB",
    "MAX",
    "TMEM127",
    "RYR2",
    "PKP2",
    "DSP",
    "DSC2",
    "TMEM43",
    "DSG2",
    "KCNJ2",
    "LMNA",
    "EMD",
    "MYBPC3",
    "MYH7",
    "TNNT2",
    "TNNI3",
    "TPM1",
    "MYL3",
    "ACTC1",
    "PRKAG2",
    "GLA",
    "MYL2",
    "LAMP2",
    "PTPN11",
    "RAF1",
    "BRAF",
    "HRAS",
    "KRAS",
    "NRAS",
    "SOS1",
    "OTC",
]


def _scrape_acmg_genes_from_ncbi(url: str) -> list[str]:
    """
    Scrape ACMG genes from NCBI ClinVar website.

    Args:
        url: URL to the NCBI ACMG page

    Returns:
        List of cleaned gene symbols

    Raises:
        ValueError: If the expected table structure is not found
        requests.RequestException: If the web request fails
    """
    logger.info(f"Attempting to scrape ACMG genes from: {url}")

    # Make request with User-Agent header
    headers = {
        "User-Agent": "custom-panel/1.0 (Python scientific tool for gene panel curation)"
    }

    try:
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
    except requests.RequestException as e:
        logger.error(f"Failed to fetch ACMG page from {url}: {e}")
        raise

    # Parse HTML content
    soup = BeautifulSoup(response.content, "html.parser")

    # Find the main content div
    main_content = soup.find("div", id="maincontent")
    if not main_content:
        raise ValueError(
            "Could not find the main content div on the NCBI page. The website structure may have changed."
        )
    assert isinstance(main_content, Tag)  # Help MyPy understand the type

    # Find the first table within main content
    table_element = main_content.find("table")
    if not table_element:
        raise ValueError(
            "Could not find the ACMG gene table on the NCBI page. The website structure may have changed."
        )

    # Use robust parsing method due to malformed HTML in the NCBI table
    # The table has issues with nested <tr> tags that break standard parsers
    logger.info("Using robust parsing method to handle malformed table HTML")

    cleaned_genes = []
    seen_genes = set()

    # Find all gene links directly in the table HTML (both GTR and regular gene links)
    # This bypasses the malformed table structure
    gene_links = table_element.find_all(
        "a", href=lambda x: x and ("/gtr/genes/" in x or "/gene/" in x)
    )

    for link in gene_links:
        gene_text = link.get_text(strip=True)
        if gene_text:
            # Remove MIM numbers and other suffixes (e.g., "COL3A1 (MIM 120180)" -> "COL3A1")
            if "(" in gene_text:
                gene_symbol = gene_text.split("(")[0].strip()
            else:
                gene_symbol = gene_text.strip()

            # Only add valid gene symbols (allow alphanumeric)
            if (
                gene_symbol
                and gene_symbol.replace("_", "").replace("-", "").isalnum()
                and gene_symbol not in seen_genes
            ):
                cleaned_genes.append(gene_symbol)
                seen_genes.add(gene_symbol)
                logger.debug(f"Extracted gene: {gene_symbol}")

    if not cleaned_genes:
        # Fallback: try pandas parsing method
        logger.warning("No genes found with robust method, trying pandas fallback")
        try:
            table_html = StringIO(str(table_element))
            tables = pd.read_html(table_html)
            if tables:
                df = tables[0]
                expected_column = "Gene via GTR"
                if expected_column in df.columns:
                    gene_series = df[expected_column].dropna()
                    for cell_content in gene_series:
                        cell_soup = BeautifulSoup(str(cell_content), "html.parser")
                        gene_link = cell_soup.find(
                            "a",
                            href=lambda x: x and ("/gtr/genes/" in x or "/gene/" in x),
                        )
                        if gene_link:
                            gene_symbol = gene_link.get_text(strip=True)
                            if "(" in gene_symbol:
                                gene_symbol = gene_symbol.split("(")[0].strip()
                            if gene_symbol and gene_symbol not in seen_genes:
                                cleaned_genes.append(gene_symbol)
                                seen_genes.add(gene_symbol)
        except Exception as e:
            logger.warning(f"Pandas fallback also failed: {e}")

    logger.info(f"Successfully scraped {len(cleaned_genes)} ACMG genes from NCBI")
    logger.debug(f"Extracted genes: {sorted(cleaned_genes)}")
    return cleaned_genes


def fetch_acmg_incidental_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch ACMG incidental findings gene data.

    Prioritizes live scraping from NCBI, with fallbacks to file and default list.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with ACMG incidental findings data
    """
    acmg_config = config.get("data_sources", {}).get("ACMG_Incidental_Findings", {})

    if not acmg_config.get("enabled", True):
        logger.info("ACMG incidental findings data source is disabled")
        return pd.DataFrame()

    evidence_score = acmg_config.get("evidence_score", 1.5)
    url = acmg_config.get("url")
    file_path = acmg_config.get("file_path")
    retrieval_date = datetime.now().strftime("%Y-%m-%d")

    genes = []
    source_details = []

    # Priority 1: Try to scrape from live URL
    if url:
        try:
            genes = _scrape_acmg_genes_from_ncbi(url)
            if genes:
                logger.info(
                    f"Successfully fetched {len(genes)} ACMG genes from live URL: {url}"
                )
                source_details = [f"URL:{url}|Date:{retrieval_date}"] * len(genes)
        except Exception as e:  # Catch all exceptions from scraper
            logger.warning(f"Failed to scrape ACMG genes from {url}: {e}")
            logger.warning("Attempting to use fallback data sources")

    # Priority 2: Fall back to file if scraping failed
    if not genes and file_path:
        try:
            if Path(file_path).exists():
                file_path_obj = Path(file_path)
                if file_path_obj.suffix.lower() == ".csv":
                    # Handle CSV files
                    df_fallback = pd.read_csv(file_path)
                    # Try common column names for gene symbols
                    gene_columns = [
                        "Gene",
                        "Gene Symbol",
                        "gene_symbol",
                        "symbol",
                        "approved_symbol",
                    ]
                    gene_column = None
                    for col in gene_columns:
                        if col in df_fallback.columns:
                            gene_column = col
                            break

                    if gene_column:
                        genes = (
                            df_fallback[gene_column]
                            .dropna()
                            .astype(str)
                            .str.strip()
                            .tolist()
                        )
                        genes = [
                            gene for gene in genes if gene and not gene.startswith("#")
                        ]
                        logger.info(
                            f"Loaded {len(genes)} ACMG genes from CSV fallback file: {file_path}"
                        )
                        source_details = [f"File:{file_path_obj.name}"] * len(genes)
                else:
                    # Handle text files (one gene per line)
                    with open(file_path, encoding="utf-8") as f:
                        lines = f.readlines()

                    genes = []
                    for line in lines:
                        line = line.strip()
                        if (
                            line
                            and not line.startswith("#")
                            and not line.startswith("//")
                        ):
                            genes.append(line)

                    if genes:
                        logger.info(
                            f"Loaded {len(genes)} ACMG genes from text fallback file: {file_path}"
                        )
                        source_details = [f"File:{file_path_obj.name}"] * len(genes)
        except Exception as e:
            logger.warning(
                f"Could not load ACMG genes from fallback file {file_path}: {e}"
            )

    # Priority 3: Use default list if all else fails
    if not genes:
        genes = DEFAULT_ACMG_GENES.copy()
        logger.warning("Using default ACMG gene list as final fallback")
        logger.info(f"Using default ACMG gene list with {len(genes)} genes")
        source_details = ["ACMG_SF_v3.2_default"] * len(genes)

    if not genes:
        logger.error("No ACMG genes available from any source")
        return pd.DataFrame()

    # Create standardized dataframe
    evidence_scores = [evidence_score] * len(genes)

    df = create_standard_dataframe(
        genes=genes,
        source_name="ACMG_Incidental_Findings",
        evidence_scores=evidence_scores,
        source_details=source_details,
        gene_names_reported=genes,
    )

    logger.info(f"Created ACMG incidental findings dataset with {len(df)} genes")
    return df


def get_acmg_gene_categories() -> dict[str, list[str]]:
    """
    Get ACMG genes organized by clinical category.

    Returns:
        Dictionary mapping categories to gene lists
    """
    categories = {
        "Hereditary_Cancer": [
            "BRCA1",
            "BRCA2",
            "TP53",
            "STK11",
            "MLH1",
            "MSH2",
            "MSH6",
            "PMS2",
            "EPCAM",
            "VHL",
            "MEN1",
            "RET",
            "PTEN",
            "TSC1",
            "TSC2",
            "WT1",
            "NF2",
            "BMPR1A",
            "SMAD4",
            "MUTYH",
        ],
        "Cardiovascular": [
            "COL3A1",
            "FBN1",
            "TGFBR1",
            "TGFBR2",
            "SMAD3",
            "ACTA2",
            "MYLK",
            "MYH11",
            "KCNQ1",
            "KCNH2",
            "SCN5A",
            "LDLR",
            "APOB",
            "PCSK9",
            "RYR2",
            "PKP2",
            "DSP",
            "DSC2",
            "TMEM43",
            "DSG2",
            "KCNJ2",
            "LMNA",
            "EMD",
            "MYBPC3",
            "MYH7",
            "TNNT2",
            "TNNI3",
            "TPM1",
            "MYL3",
            "ACTC1",
            "PRKAG2",
            "MYL2",
            "LAMP2",
        ],
        "Metabolic": ["RYR1", "CACNA1S", "ATP7B", "GLA", "OTC"],
        "Endocrine": ["SDHD", "SDHAF2", "SDHC", "SDHB", "MAX", "TMEM127"],
        "RASopathy": ["PTPN11", "RAF1", "BRAF", "HRAS", "KRAS", "NRAS", "SOS1"],
    }

    return categories


def validate_acmg_config(config: dict[str, Any]) -> list[str]:
    """
    Validate ACMG incidental findings configuration.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation errors
    """
    errors = []
    acmg_config = config.get("data_sources", {}).get("ACMG_Incidental_Findings", {})

    # Check evidence score
    evidence_score = acmg_config.get("evidence_score", 1.5)
    if not isinstance(evidence_score, int | float) or evidence_score < 0:
        errors.append("ACMG evidence score must be a non-negative number")

    # Check file path if specified
    file_path = acmg_config.get("file_path")
    if file_path and not Path(file_path).exists():
        errors.append(f"ACMG file not found: {file_path}")

    return errors


def get_acmg_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary of ACMG incidental findings configuration.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    acmg_config = config.get("data_sources", {}).get("ACMG_Incidental_Findings", {})

    summary = {
        "enabled": acmg_config.get("enabled", True),
        "evidence_score": acmg_config.get("evidence_score", 1.5),
        "file_path": acmg_config.get("file_path"),
        "file_exists": False,
        "gene_count": 0,
        "categories": get_acmg_gene_categories(),
        "validation_errors": validate_acmg_config(config),
    }

    # Check file and count genes
    file_path = acmg_config.get("file_path")
    if file_path:
        summary["file_exists"] = Path(file_path).exists()

    # Always show default count in summary - actual count depends on runtime scraping
    summary["gene_count"] = len(DEFAULT_ACMG_GENES)
    summary["using_default"] = True  # This will be updated at runtime
    summary["url"] = acmg_config.get("url")

    return summary
