"""
HPO/OMIM neoplasm data source extractor.

This module identifies genes associated with neoplasia by processing HPO
(Human Phenotype Ontology) and OMIM data through API interactions and file parsing.
"""

import logging
import re
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.config_manager import ConfigManager
from ..core.hpo_client import HPOClient
from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)

# HPO root term for neoplasm
NEOPLASM_ROOT_TERM = "HP:0002664"

# URLs for data files
PHENOTYPE_HPOA_URL = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"


def fetch_hpo_neoplasm_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch neoplasm-associated genes from HPO/OMIM data sources.

    This follows the R implementation logic:
    1. Get all HPO descendants of the neoplasm term (HP:0002664)
    2. Download and parse phenotype.hpoa file
    3. Filter phenotype.hpoa for neoplasm HPO terms
    4. Join with OMIM genemap2 data to get associated genes
    5. Return standardized gene list

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with HPO neoplasm data
    """
    config_manager = ConfigManager(config)
    hpo_config = config_manager.get_source_config("HPO_Neoplasm")

    if not hpo_config.get("enabled", True):
        logger.info("HPO neoplasm data source is disabled")
        return pd.DataFrame()

    # Initialize HPO client
    client = HPOClient()

    # Get cache configuration
    cache_dir = Path(hpo_config.get("cache_dir", ".cache/hpo"))
    cache_expiry_days = hpo_config.get("cache_expiry_days", 30)

    try:
        # Step 1: Get all descendants of neoplasm term
        logger.info(f"Fetching descendants of {NEOPLASM_ROOT_TERM}")
        neoplasm_terms = client.get_descendant_terms(
            NEOPLASM_ROOT_TERM, max_depth=20, include_self=True
        )
        logger.info(f"Found {len(neoplasm_terms)} neoplasm-related HPO terms")

        # Step 2: Download/load phenotype.hpoa file
        phenotype_hpoa_path = cache_dir / "phenotype.hpoa"
        if not client.is_cache_valid(phenotype_hpoa_path, cache_expiry_days):
            try:
                client.download_file(PHENOTYPE_HPOA_URL, phenotype_hpoa_path)
            except Exception as e:
                logger.error(f"Failed to download phenotype.hpoa: {e}")
                if phenotype_hpoa_path.exists():
                    logger.warning("Using expired cache file")
                else:
                    logger.error("No phenotype.hpoa cache available")
                    return pd.DataFrame()
        else:
            logger.info("Using cached phenotype.hpoa file")

        # Parse phenotype.hpoa
        phenotype_df = client.parse_phenotype_hpoa(phenotype_hpoa_path)

        # Step 3: Filter for neoplasm HPO terms
        neoplasm_phenotypes = phenotype_df[phenotype_df["hpo_id"].isin(neoplasm_terms)]
        logger.info(f"Found {len(neoplasm_phenotypes)} neoplasm phenotype annotations")

        # Extract unique disease IDs (OMIM IDs)
        omim_diseases = neoplasm_phenotypes[
            neoplasm_phenotypes["database_id"].str.startswith("OMIM:")
        ]
        unique_omim_ids = (
            omim_diseases["database_id"].str.extract(r"OMIM:(\d+)")[0].unique()
        )
        logger.info(f"Found {len(unique_omim_ids)} unique OMIM disease IDs")

        # Step 4: Load OMIM genemap2 data (from URL or local file)
        omim_genemap2_url = hpo_config.get("omim_genemap2_url")
        omim_genemap2_path = hpo_config.get("omim_genemap2_path")

        # Debug logging to trace configuration
        logger.info(f"DEBUG: HPO config - omim_genemap2_url: {omim_genemap2_url}")
        logger.info(f"DEBUG: HPO config - omim_genemap2_path: {omim_genemap2_path}")

        if omim_genemap2_url:
            # Download from URL
            logger.info(f"DEBUG: Using OMIM genemap2 URL: {omim_genemap2_url}")
            cache_dir = Path(hpo_config.get("cache_dir", ".cache/hpo"))
            cache_dir.mkdir(parents=True, exist_ok=True)

            genemap2_cache_path = cache_dir / "genemap2.txt"
            cache_expiry_days = hpo_config.get("cache_expiry_days", 30)
            logger.info(f"DEBUG: Cache path will be: {genemap2_cache_path}")

            if not client.is_cache_valid(genemap2_cache_path, cache_expiry_days):
                logger.info("Downloading OMIM genemap2 from URL")
                try:
                    client.download_file(omim_genemap2_url, genemap2_cache_path)
                    logger.info(
                        f"DEBUG: Successfully downloaded to cache: {genemap2_cache_path}"
                    )
                except Exception as e:
                    logger.error(f"Failed to download OMIM genemap2: {e}")
                    if genemap2_cache_path.exists():
                        logger.warning("Using expired cache file")
                    else:
                        logger.error("No OMIM genemap2 cache available")
                        return pd.DataFrame()
            else:
                logger.info(
                    f"DEBUG: Using cached OMIM genemap2 file: {genemap2_cache_path}"
                )

            omim_df = client.parse_omim_genemap2(genemap2_cache_path)
            logger.info(
                f"DEBUG: Parsed {len(omim_df)} entries from cached file: {genemap2_cache_path}"
            )

        elif omim_genemap2_path:
            # Use local file
            logger.info(f"DEBUG: Using OMIM genemap2 local path: {omim_genemap2_path}")
            omim_path = Path(omim_genemap2_path)
            if not omim_path.exists():
                logger.error(f"OMIM genemap2 file not found: {omim_path}")
                return pd.DataFrame()

            omim_df = client.parse_omim_genemap2(omim_path)
            logger.info(
                f"DEBUG: Parsed {len(omim_df)} entries from local file: {omim_path}"
            )

        else:
            logger.error("DEBUG: No OMIM genemap2 URL or file path configured")
            logger.error("Please set 'omim_genemap2_url' in your config.local.yml file")
            logger.error("Get your access token from: https://omim.org/downloads")
            return pd.DataFrame()

        # Step 5: Extract genes associated with neoplasm OMIM IDs
        all_genes = _extract_genes_from_omim(
            omim_df, unique_omim_ids.tolist(), neoplasm_phenotypes, hpo_config
        )

        if not all_genes:
            logger.warning("No genes found associated with neoplasm phenotypes")
            return pd.DataFrame()

        # Convert to standardized DataFrame
        genes_list = []
        evidence_scores = []
        source_details = []

        for gene_symbol, gene_data in all_genes.items():
            genes_list.append(gene_symbol)

            # Calculate evidence score
            evidence_score = _calculate_evidence_score(gene_data, hpo_config)
            evidence_scores.append(evidence_score)

            # Create source details
            details = _create_source_details(gene_data)
            source_details.append(details)

        standardized_df = create_standard_dataframe(
            genes=genes_list,
            source_name="HPO_Neoplasm",
            evidence_scores=evidence_scores,
            source_details=source_details,
            gene_names_reported=genes_list,
        )

        logger.info(
            f"Successfully processed {len(standardized_df)} HPO/OMIM neoplasm genes"
        )
        return standardized_df

    except Exception as e:
        logger.error(f"Error processing HPO neoplasm data: {e}")
        return pd.DataFrame()


def _extract_genes_from_omim(
    omim_df: pd.DataFrame,
    omim_ids: list[str],
    phenotype_df: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """
    Extract genes associated with specific OMIM IDs.

    Args:
        omim_df: OMIM genemap2 DataFrame
        omim_ids: List of OMIM IDs to search for
        phenotype_df: Phenotype annotations DataFrame
        config: Configuration dictionary

    Returns:
        Dictionary mapping gene symbols to gene data
    """
    all_genes: dict[str, dict[str, Any]] = {}

    # Create mapping of OMIM ID to phenotype data
    omim_to_phenotypes = {}
    for omim_id in omim_ids:
        full_omim_id = f"OMIM:{omim_id}"
        related_phenotypes = phenotype_df[phenotype_df["database_id"] == full_omim_id]
        if not related_phenotypes.empty:
            omim_to_phenotypes[omim_id] = related_phenotypes

    logger.info(f"Processing {len(omim_to_phenotypes)} OMIM IDs with phenotypes")

    # Process each row in OMIM genemap2
    for _, row in omim_df.iterrows():
        # Check if this row has phenotypes
        phenotypes_str = str(row.get("phenotypes", ""))
        if pd.isna(phenotypes_str) or not phenotypes_str:
            continue

        # Extract OMIM IDs from phenotypes column using regex
        # Format: "Phenotype name, MIM_NUMBER (inheritance), ..."
        omim_pattern = re.compile(r"\b(\d{6})\b")
        row_omim_ids = set(omim_pattern.findall(phenotypes_str))

        # Check if any of our target OMIM IDs match
        matching_omim_ids = row_omim_ids.intersection(set(omim_ids))
        if not matching_omim_ids:
            continue

        # Extract gene symbol
        gene_symbol = row.get("approved_symbol", "")
        if pd.isna(gene_symbol) or not gene_symbol:
            gene_symbol = row.get("gene_symbols", "")
            if pd.isna(gene_symbol) or not gene_symbol:
                continue

        # Handle multiple gene symbols
        gene_symbols = [s.strip() for s in str(gene_symbol).split(",")]

        for symbol in gene_symbols:
            if not symbol:
                continue

            symbol = symbol.upper()

            if symbol not in all_genes:
                all_genes[symbol] = {
                    "gene_symbol": symbol,
                    "entrez_id": row.get("entrez_gene_id"),
                    "mim_number": row.get("mim_number"),
                    "hpo_terms": set(),
                    "diseases": set(),
                    "omim_ids": set(),
                    "phenotypes": [],
                    "evidence_sources": {"HPO", "OMIM"},
                }

            # Add data from matching OMIM IDs
            for omim_id in matching_omim_ids:
                all_genes[symbol]["omim_ids"].add(omim_id)

                # Add associated HPO terms and disease names
                if omim_id in omim_to_phenotypes:
                    phenotypes = omim_to_phenotypes[omim_id]
                    for _, pheno_row in phenotypes.iterrows():
                        all_genes[symbol]["hpo_terms"].add(pheno_row["hpo_id"])
                        all_genes[symbol]["diseases"].add(pheno_row["disease_name"])
                        all_genes[symbol]["phenotypes"].append(
                            {
                                "hpo_id": pheno_row["hpo_id"],
                                "disease_name": pheno_row["disease_name"],
                                "evidence": pheno_row.get("evidence", ""),
                                "frequency": pheno_row.get("frequency", ""),
                            }
                        )

    logger.info(f"Found {len(all_genes)} genes associated with neoplasm phenotypes")
    return all_genes


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
    base_score = config.get("base_evidence_score", 0.7)

    # Bonus for multiple HPO terms
    hpo_count = len(gene_data.get("hpo_terms", set()))
    hpo_bonus = min(hpo_count * 0.02, 0.15)

    # Bonus for multiple diseases
    disease_count = len(gene_data.get("diseases", set()))
    disease_bonus = min(disease_count * 0.02, 0.1)

    # Bonus for having Entrez ID
    entrez_bonus = 0.05 if gene_data.get("entrez_id") else 0.0

    final_score = base_score + hpo_bonus + disease_bonus + entrez_bonus
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

    # HPO terms count
    hpo_count = len(gene_data.get("hpo_terms", set()))
    if hpo_count > 0:
        details.append(f"HPO_terms:{hpo_count}")

    # Diseases count
    disease_count = len(gene_data.get("diseases", set()))
    if disease_count > 0:
        details.append(f"Diseases:{disease_count}")

    # OMIM IDs
    omim_ids = gene_data.get("omim_ids", set())
    if omim_ids:
        details.append(f"OMIM_IDs:{len(omim_ids)}")

    # MIM number
    mim_number = gene_data.get("mim_number")
    if mim_number:
        details.append(f"MIM:{mim_number}")

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
    errors: list[str] = []
    config_manager = ConfigManager(config)
    hpo_config = config_manager.get_source_config("HPO_Neoplasm")

    if not isinstance(hpo_config, dict):
        errors.append("hpo_neoplasm config must be a dictionary")
        return errors

    # Check if enabled
    if not hpo_config.get("enabled", True):
        return errors  # Skip further validation if disabled

    # Validate OMIM genemap2 file path (required)
    omim_file = hpo_config.get("omim_genemap2_path")
    if not omim_file:
        errors.append("omim_genemap2_path is required for HPO neoplasm data source")
    elif not isinstance(omim_file, str):
        errors.append("omim_genemap2_path must be a string")
    elif not Path(omim_file).exists():
        errors.append(f"OMIM genemap2 file does not exist: {omim_file}")

    # Validate cache settings
    cache_expiry = hpo_config.get("cache_expiry_days")
    if cache_expiry is not None:
        if not isinstance(cache_expiry, int) or cache_expiry < 1:
            errors.append("cache_expiry_days must be a positive integer")

    # Validate evidence score
    base_score = hpo_config.get("base_evidence_score")
    if base_score is not None:
        if not isinstance(base_score, int | float) or not 0 <= base_score <= 1:
            errors.append("base_evidence_score must be between 0.0 and 1.0")

    return errors


def get_hpo_neoplasm_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary of HPO neoplasm configuration.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    config_manager = ConfigManager(config)
    hpo_config = config_manager.get_source_config("HPO_Neoplasm")

    summary = {
        "enabled": hpo_config.get("enabled", True),
        "omim_genemap2_configured": bool(hpo_config.get("omim_genemap2_path")),
        "cache_dir": hpo_config.get("cache_dir", ".cache/hpo"),
        "cache_expiry_days": hpo_config.get("cache_expiry_days", 30),
        "base_evidence_score": hpo_config.get("base_evidence_score", 0.7),
        "validation_errors": validate_hpo_neoplasm_config(config),
    }

    # Check cache status
    cache_dir = Path(summary["cache_dir"])
    phenotype_hpoa_path = cache_dir / "phenotype.hpoa"

    if phenotype_hpoa_path.exists():
        summary["phenotype_hpoa_cached"] = True
        summary["phenotype_hpoa_cache_valid"] = HPOClient.is_cache_valid(
            phenotype_hpoa_path, summary["cache_expiry_days"]
        )
    else:
        summary["phenotype_hpoa_cached"] = False
        summary["phenotype_hpoa_cache_valid"] = False

    # Add download URL info
    summary["phenotype_hpoa_url"] = PHENOTYPE_HPOA_URL
    summary["neoplasm_root_term"] = NEOPLASM_ROOT_TERM

    return summary
