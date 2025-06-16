"""
ACMG Incidental Findings data source extractor.

This module processes ACMG recommendations for reporting incidental findings
in clinical exome and genome sequencing.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

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


def fetch_acmg_incidental_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch ACMG incidental findings gene data.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with ACMG incidental findings data
    """
    acmg_config = config.get("data_sources", {}).get("acmg_incidental", {})

    if not acmg_config.get("enabled", True):
        logger.info("ACMG incidental findings data source is disabled")
        return pd.DataFrame()

    evidence_score = acmg_config.get("evidence_score", 1.5)
    file_path = acmg_config.get("file_path")

    # Try to load from file first, then fall back to default list
    genes = []
    source_details = []

    if file_path:
        genes = load_acmg_genes_from_file(file_path)
        if genes:
            logger.info(f"Loaded {len(genes)} ACMG genes from file: {file_path}")
            source_details = [f"File:{Path(file_path).name}"] * len(genes)
        else:
            logger.warning(
                f"Could not load ACMG genes from file {file_path}, using default list"
            )

    # Use default list if no file or file loading failed
    if not genes:
        genes = DEFAULT_ACMG_GENES.copy()
        logger.info(f"Using default ACMG gene list with {len(genes)} genes")
        source_details = ["ACMG_SF_v3.2_default"] * len(genes)

    if not genes:
        logger.warning("No ACMG genes available")
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


def load_acmg_genes_from_file(file_path: str | Path) -> list[str]:
    """
    Load ACMG genes from a file.

    Args:
        file_path: Path to the file containing ACMG genes

    Returns:
        List of gene symbols
    """
    file_path = Path(file_path)

    if not file_path.exists():
        logger.error(f"ACMG file not found: {file_path}")
        return []

    try:
        if file_path.suffix.lower() in [".xlsx", ".xls"]:
            return load_acmg_genes_from_excel(file_path)
        elif file_path.suffix.lower() == ".csv":
            return load_acmg_genes_from_csv(file_path)
        elif file_path.suffix.lower() in [".txt", ".tsv"]:
            return load_acmg_genes_from_text(file_path)
        else:
            logger.error(f"Unsupported file format for ACMG genes: {file_path.suffix}")
            return []
    except Exception as e:
        logger.error(f"Error loading ACMG genes from {file_path}: {e}")
        return []


def load_acmg_genes_from_excel(file_path: Path) -> list[str]:
    """Load ACMG genes from Excel file."""
    try:
        # Try different possible column names and sheets
        possible_columns = [
            "Gene",
            "Gene Symbol",
            "Symbol",
            "Gene_Symbol",
            "GENE",
            "gene",
        ]

        # Try first sheet
        df = pd.read_excel(file_path, engine="openpyxl")

        # Find gene column
        gene_column = None
        for col in possible_columns:
            if col in df.columns:
                gene_column = col
                break

        if not gene_column:
            # Try the first column if no standard gene column found
            if len(df.columns) > 0:
                gene_column = df.columns[0]
                logger.warning(f"Using first column '{gene_column}' as gene column")
            else:
                logger.error("No columns found in Excel file")
                return []

        genes = df[gene_column].dropna().astype(str).str.strip().tolist()
        genes = [gene for gene in genes if gene]  # Remove empty strings

        return list(dict.fromkeys(genes))  # Remove duplicates while preserving order
    except Exception as e:
        logger.error(f"Error reading Excel file {file_path}: {e}")
        return []


def load_acmg_genes_from_csv(file_path: Path) -> list[str]:
    """Load ACMG genes from CSV file."""
    try:
        # Try different encodings
        for encoding in ["utf-8", "latin1", "cp1252"]:
            try:
                df = pd.read_csv(file_path, encoding=encoding)
                break
            except UnicodeDecodeError:
                continue
        else:
            logger.error(f"Could not read CSV file {file_path} with any encoding")
            return []

        # Find gene column (similar to Excel method)
        possible_columns = [
            "Gene",
            "Gene Symbol",
            "Symbol",
            "Gene_Symbol",
            "GENE",
            "gene",
        ]

        gene_column = None
        for col in possible_columns:
            if col in df.columns:
                gene_column = col
                break

        if not gene_column:
            if len(df.columns) > 0:
                gene_column = df.columns[0]
                logger.warning(f"Using first column '{gene_column}' as gene column")
            else:
                logger.error("No columns found in CSV file")
                return []

        genes = df[gene_column].dropna().astype(str).str.strip().tolist()
        genes = [gene for gene in genes if gene]

        return list(dict.fromkeys(genes))
    except Exception as e:
        logger.error(f"Error reading CSV file {file_path}: {e}")
        return []


def load_acmg_genes_from_text(file_path: Path) -> list[str]:
    """Load ACMG genes from text file (one gene per line)."""
    try:
        with open(file_path, encoding="utf-8") as f:
            genes = [line.strip() for line in f if line.strip()]

        # Handle tab-separated files
        if file_path.suffix.lower() == ".tsv":
            # Assume genes are in the first column
            processed_genes = []
            for line in genes:
                if "\t" in line:
                    processed_genes.append(line.split("\t")[0].strip())
                else:
                    processed_genes.append(line)
            genes = processed_genes

        genes = [
            gene for gene in genes if gene and not gene.startswith("#")
        ]  # Remove comments
        return list(dict.fromkeys(genes))
    except Exception as e:
        logger.error(f"Error reading text file {file_path}: {e}")
        return []


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
    acmg_config = config.get("data_sources", {}).get("acmg_incidental", {})

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
    acmg_config = config.get("data_sources", {}).get("acmg_incidental", {})

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
        if summary["file_exists"]:
            genes = load_acmg_genes_from_file(file_path)
            summary["gene_count"] = len(genes)

    if summary["gene_count"] == 0:
        summary["gene_count"] = len(DEFAULT_ACMG_GENES)
        summary["using_default"] = True
    else:
        summary["using_default"] = False

    return summary
