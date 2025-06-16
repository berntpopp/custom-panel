"""
HPO/OMIM neoplasm data source extractor.

This module identifies genes associated with neoplasia by processing HPO
(Human Phenotype Ontology) and OMIM data through API interactions and file parsing.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.hpo_client import HPOClient
from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)


def fetch_hpo_neoplasm_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch neoplasm-associated genes from HPO/OMIM data sources.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with HPO neoplasm data
    """
    hpo_config = config.get("data_sources", {}).get("hpo_neoplasm", {})

    if not hpo_config.get("enabled", True):
        logger.info("HPO neoplasm data source is disabled")
        return pd.DataFrame()

    # Initialize HPO client
    api_config = config.get("apis", {}).get("hpo", {})
    client = HPOClient(
        timeout=api_config.get("timeout", 30),
        max_retries=api_config.get("max_retries", 3),
        retry_delay=api_config.get("retry_delay", 1.0),
    )

    all_genes = {}

    # Method 1: Search for neoplasm terms and collect genes
    if hpo_config.get("use_neoplasm_search", True):
        neoplasm_genes = _get_genes_from_neoplasm_search(client, hpo_config)
        _merge_gene_data(all_genes, neoplasm_genes, "neoplasm_search")

    # Method 2: Use specific HPO terms if provided
    specific_terms = hpo_config.get("specific_hpo_terms", [])
    if specific_terms:
        specific_genes = _get_genes_from_specific_terms(
            client, specific_terms, hpo_config
        )
        _merge_gene_data(all_genes, specific_genes, "specific_terms")

    # Method 3: Load from OMIM file if provided
    omim_file = hpo_config.get("omim_file_path")
    if omim_file:
        omim_genes = _get_genes_from_omim_file(omim_file, hpo_config)
        _merge_gene_data(all_genes, omim_genes, "omim_file")

    if not all_genes:
        logger.warning("No HPO neoplasm genes found")
        return pd.DataFrame()

    # Convert to standardized DataFrame
    genes_list = []
    evidence_scores = []
    source_details = []

    for gene_symbol, gene_data in all_genes.items():
        genes_list.append(gene_symbol)

        # Calculate evidence score based on number of sources and confidence
        evidence_score = _calculate_evidence_score(gene_data, hpo_config)
        evidence_scores.append(evidence_score)

        # Create source details
        details = _create_source_details(gene_data)
        source_details.append(details)

    standardized_df = create_standard_dataframe(
        genes=genes_list,
        source_name="HPO_OMIM_Neoplasm",
        evidence_scores=evidence_scores,
        source_details=source_details,
        gene_names_reported=genes_list,
    )

    logger.info(f"Successfully processed {len(standardized_df)} HPO neoplasm genes")
    return standardized_df


def _get_genes_from_neoplasm_search(
    client: HPOClient, config: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    """
    Get genes by searching for neoplasm-related HPO terms.

    Args:
        client: HPO client instance
        config: HPO configuration

    Returns:
        Dictionary of gene data
    """
    max_hierarchy_depth = config.get("max_hierarchy_depth", 5)

    # Find neoplasm terms
    neoplasm_terms = client.find_neoplasm_terms()
    logger.info(f"Found {len(neoplasm_terms)} neoplasm HPO terms")

    all_genes = {}

    # Get genes for each neoplasm term and its descendants
    for term in neoplasm_terms:
        term_id = term.get("id")
        if not term_id:
            continue

        logger.debug(f"Processing neoplasm term: {term_id} - {term.get('name')}")

        # Get genes for this term and its hierarchy
        hierarchy_genes = client.get_genes_for_phenotype_hierarchy(
            term_id, max_depth=max_hierarchy_depth
        )

        # Merge with existing genes
        for gene_symbol, gene_data in hierarchy_genes.items():
            if gene_symbol not in all_genes:
                all_genes[gene_symbol] = {
                    "gene_symbol": gene_symbol,
                    "entrez_id": gene_data.get("entrez_id"),
                    "hpo_terms": [],
                    "diseases": set(),
                    "evidence_sources": set(),
                    "source_methods": set(),
                }

            # Merge data
            all_genes[gene_symbol]["hpo_terms"].extend(gene_data.get("hpo_terms", []))
            all_genes[gene_symbol]["diseases"].update(gene_data.get("diseases", []))
            all_genes[gene_symbol]["evidence_sources"].update(
                gene_data.get("evidence_sources", [])
            )
            all_genes[gene_symbol]["source_methods"].add("neoplasm_search")

    logger.info(f"Found {len(all_genes)} genes from neoplasm search")
    return all_genes


def _get_genes_from_specific_terms(
    client: HPOClient, hpo_terms: list[str], config: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    """
    Get genes from specific HPO terms.

    Args:
        client: HPO client instance
        hpo_terms: List of specific HPO term IDs
        config: HPO configuration

    Returns:
        Dictionary of gene data
    """
    max_hierarchy_depth = config.get("max_hierarchy_depth", 5)
    all_genes = {}

    for term_id in hpo_terms:
        logger.debug(f"Processing specific HPO term: {term_id}")

        # Get genes for this term and its hierarchy
        hierarchy_genes = client.get_genes_for_phenotype_hierarchy(
            term_id, max_depth=max_hierarchy_depth
        )

        # Merge with existing genes
        for gene_symbol, gene_data in hierarchy_genes.items():
            if gene_symbol not in all_genes:
                all_genes[gene_symbol] = {
                    "gene_symbol": gene_symbol,
                    "entrez_id": gene_data.get("entrez_id"),
                    "hpo_terms": [],
                    "diseases": set(),
                    "evidence_sources": set(),
                    "source_methods": set(),
                }

            # Merge data
            all_genes[gene_symbol]["hpo_terms"].extend(gene_data.get("hpo_terms", []))
            all_genes[gene_symbol]["diseases"].update(gene_data.get("diseases", []))
            all_genes[gene_symbol]["evidence_sources"].update(
                gene_data.get("evidence_sources", [])
            )
            all_genes[gene_symbol]["source_methods"].add("specific_terms")

    logger.info(f"Found {len(all_genes)} genes from specific HPO terms")
    return all_genes


def _get_genes_from_omim_file(
    file_path: str, config: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    """
    Get genes from OMIM file.

    Args:
        file_path: Path to OMIM file
        config: HPO configuration

    Returns:
        Dictionary of gene data
    """
    file_path = Path(file_path)
    if not file_path.exists():
        logger.error(f"OMIM file not found: {file_path}")
        return {}

    try:
        # Determine file format and read
        if file_path.suffix.lower() == ".csv":
            df = pd.read_csv(file_path)
        elif file_path.suffix.lower() in [".xlsx", ".xls"]:
            df = pd.read_excel(file_path)
        elif file_path.suffix.lower() == ".txt":
            # Assume tab-separated for .txt files
            df = pd.read_csv(file_path, sep="\t")
        else:
            logger.error(f"Unsupported OMIM file format: {file_path.suffix}")
            return {}

        if df.empty:
            logger.warning(f"OMIM file is empty: {file_path}")
            return {}

        # Extract gene information
        gene_column = config.get("omim_gene_column", "Gene Symbol")
        disease_column = config.get("omim_disease_column", "Disease")
        omim_id_column = config.get("omim_id_column", "OMIM ID")

        if gene_column not in df.columns:
            logger.error(f"Gene column '{gene_column}' not found in OMIM file")
            return {}

        all_genes = {}

        for _, row in df.iterrows():
            gene_symbol = row.get(gene_column)
            if pd.isna(gene_symbol) or not gene_symbol:
                continue

            gene_symbol = str(gene_symbol).strip().upper()

            if gene_symbol not in all_genes:
                all_genes[gene_symbol] = {
                    "gene_symbol": gene_symbol,
                    "entrez_id": None,
                    "hpo_terms": [],
                    "diseases": set(),
                    "evidence_sources": set(),
                    "source_methods": set(),
                    "omim_data": [],
                }

            # Add OMIM data
            omim_entry = {
                "disease": row.get(disease_column, ""),
                "omim_id": row.get(omim_id_column, ""),
            }
            all_genes[gene_symbol]["omim_data"].append(omim_entry)

            if row.get(disease_column):
                all_genes[gene_symbol]["diseases"].add(str(row.get(disease_column)))

            all_genes[gene_symbol]["evidence_sources"].add("OMIM")
            all_genes[gene_symbol]["source_methods"].add("omim_file")

        logger.info(f"Found {len(all_genes)} genes from OMIM file")
        return all_genes

    except Exception as e:
        logger.error(f"Error reading OMIM file {file_path}: {e}")
        return {}


def _merge_gene_data(
    target: dict[str, dict[str, Any]], source: dict[str, dict[str, Any]], method: str
) -> None:
    """
    Merge gene data from source into target dictionary.

    Args:
        target: Target gene dictionary to merge into
        source: Source gene dictionary to merge from
        method: Method name for tracking data sources
    """
    for gene_symbol, gene_data in source.items():
        if gene_symbol not in target:
            target[gene_symbol] = {
                "gene_symbol": gene_symbol,
                "entrez_id": gene_data.get("entrez_id"),
                "hpo_terms": [],
                "diseases": set(),
                "evidence_sources": set(),
                "source_methods": set(),
            }
            # Copy OMIM data if present
            if "omim_data" in gene_data:
                target[gene_symbol]["omim_data"] = gene_data["omim_data"]

        # Merge data
        target[gene_symbol]["hpo_terms"].extend(gene_data.get("hpo_terms", []))
        target[gene_symbol]["diseases"].update(gene_data.get("diseases", set()))
        target[gene_symbol]["evidence_sources"].update(
            gene_data.get("evidence_sources", set())
        )
        target[gene_symbol]["source_methods"].add(method)

        # Update entrez_id if not set
        if not target[gene_symbol]["entrez_id"] and gene_data.get("entrez_id"):
            target[gene_symbol]["entrez_id"] = gene_data.get("entrez_id")


def _calculate_evidence_score(
    gene_data: dict[str, Any], config: dict[str, Any]
) -> float:
    """
    Calculate evidence score for a gene based on available data.

    Args:
        gene_data: Gene data dictionary
        config: HPO configuration

    Returns:
        Evidence score between 0.0 and 1.0
    """
    base_score = config.get("base_evidence_score", 0.5)

    # Bonus for multiple data sources
    source_count = len(gene_data.get("source_methods", set()))
    source_bonus = min(source_count * 0.1, 0.3)

    # Bonus for multiple HPO terms
    hpo_count = len(gene_data.get("hpo_terms", []))
    hpo_bonus = min(hpo_count * 0.05, 0.2)

    # Bonus for multiple diseases
    disease_count = len(gene_data.get("diseases", set()))
    disease_bonus = min(disease_count * 0.02, 0.1)

    # Bonus for having Entrez ID
    entrez_bonus = 0.05 if gene_data.get("entrez_id") else 0.0

    final_score = base_score + source_bonus + hpo_bonus + disease_bonus + entrez_bonus
    return min(final_score, 1.0)


def _create_source_details(gene_data: dict[str, Any]) -> str:
    """
    Create source details string for a gene.

    Args:
        gene_data: Gene data dictionary

    Returns:
        Source details string
    """
    details = []

    # Source methods
    methods = gene_data.get("source_methods", set())
    if methods:
        details.append(f"Methods:{','.join(sorted(methods))}")

    # HPO terms count
    hpo_count = len(gene_data.get("hpo_terms", []))
    if hpo_count > 0:
        details.append(f"HPO_terms:{hpo_count}")

    # Diseases count
    disease_count = len(gene_data.get("diseases", set()))
    if disease_count > 0:
        details.append(f"Diseases:{disease_count}")

    # Evidence sources
    evidence_sources = gene_data.get("evidence_sources", set())
    if evidence_sources:
        details.append(f"Evidence:{','.join(sorted(evidence_sources))}")

    # Entrez ID
    entrez_id = gene_data.get("entrez_id")
    if entrez_id:
        details.append(f"Entrez:{entrez_id}")

    return "|".join(details)


def validate_hpo_neoplasm_config(config: dict[str, Any]) -> list[str]:
    """
    Validate HPO neoplasm configuration.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation errors
    """
    errors = []
    hpo_config = config.get("data_sources", {}).get("hpo_neoplasm", {})

    if not isinstance(hpo_config, dict):
        errors.append("hpo_neoplasm config must be a dictionary")
        return errors

    # Validate OMIM file path if provided
    omim_file = hpo_config.get("omim_file_path")
    if omim_file:
        if not isinstance(omim_file, str):
            errors.append("omim_file_path must be a string")
        elif not Path(omim_file).exists():
            errors.append(f"OMIM file does not exist: {omim_file}")

    # Validate specific HPO terms if provided
    specific_terms = hpo_config.get("specific_hpo_terms", [])
    if specific_terms:
        if not isinstance(specific_terms, list):
            errors.append("specific_hpo_terms must be a list")
        else:
            for i, term in enumerate(specific_terms):
                if not isinstance(term, str) or not term.startswith("HP:"):
                    errors.append(
                        f"specific_hpo_terms[{i}] must be a valid HPO ID (HP:XXXXXXX)"
                    )

    # Validate numeric parameters
    numeric_params = {
        "base_evidence_score": (0.0, 1.0),
        "max_hierarchy_depth": (1, 20),
    }

    for param, (min_val, max_val) in numeric_params.items():
        value = hpo_config.get(param)
        if value is not None:
            if not isinstance(value, int | float):
                errors.append(f"{param} must be a number")
            elif not min_val <= value <= max_val:
                errors.append(f"{param} must be between {min_val} and {max_val}")

    return errors


def get_hpo_neoplasm_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary of HPO neoplasm configuration.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    hpo_config = config.get("data_sources", {}).get("hpo_neoplasm", {})

    summary = {
        "enabled": hpo_config.get("enabled", True),
        "use_neoplasm_search": hpo_config.get("use_neoplasm_search", True),
        "specific_hpo_terms_count": len(hpo_config.get("specific_hpo_terms", [])),
        "omim_file_configured": bool(hpo_config.get("omim_file_path")),
        "validation_errors": validate_hpo_neoplasm_config(config),
    }

    # Try to estimate total genes (requires API calls)
    if summary["enabled"]:
        try:
            # This is a rough estimate and would require actual API calls
            estimated_genes = 0
            if summary["use_neoplasm_search"]:
                estimated_genes += 500  # Rough estimate
            if summary["specific_hpo_terms_count"] > 0:
                estimated_genes += (
                    summary["specific_hpo_terms_count"] * 50
                )  # Rough estimate

            summary["estimated_genes"] = estimated_genes
        except Exception:
            summary["estimated_genes"] = "Unknown"

    return summary
