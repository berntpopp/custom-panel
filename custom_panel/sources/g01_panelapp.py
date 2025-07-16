"""
PanelApp data source extractor.

This module fetches gene panel data from PanelApp APIs (UK Genomics England and Australia).
"""

import logging
from typing import Any

import pandas as pd
import requests

from ..core.cache_manager import CacheManager
from ..core.config_manager import ConfigManager
from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)


class PanelAppClient:
    """Client for fetching data from PanelApp APIs with caching support."""

    def __init__(
        self,
        timeout: int = 30,
        max_retries: int = 3,
        cache_manager: CacheManager | None = None,
    ):
        """
        Initialize the PanelApp client.

        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            cache_manager: Optional cache manager for caching responses
        """
        self.timeout = timeout
        self.max_retries = max_retries
        self.cache_manager = cache_manager
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "custom-panel/0.1.0"})

    def _make_request(self, url: str) -> dict[str, Any] | None:
        """
        Make a request to PanelApp API with caching and retry logic.

        Args:
            url: API URL

        Returns:
            JSON response or None if failed
        """
        # Check cache first if available
        if self.cache_manager and self.cache_manager.enabled:
            cached_data = self.cache_manager.get("panelapp", url, "GET", None)
            if cached_data is not None:
                logger.info(f"Using cached PanelApp data from {url}")
                return cached_data

        for attempt in range(self.max_retries + 1):
            try:
                response = self.session.get(url, timeout=self.timeout)
                response.raise_for_status()
                json_response = response.json()

                # Cache the successful response
                if self.cache_manager and self.cache_manager.enabled:
                    self.cache_manager.set("panelapp", url, "GET", None, json_response)
                    logger.info("Cached PanelApp response for future use")

                return json_response
            except requests.RequestException as e:
                if attempt == self.max_retries:
                    logger.error(
                        f"Failed to fetch {url} after {self.max_retries} retries: {e}"
                    )
                    return None
                logger.warning(
                    f"Request failed (attempt {attempt + 1}/{self.max_retries + 1}): {e}"
                )
        return None

    def get_panel_list(self, base_url: str) -> list[dict[str, Any]]:
        """
        Get list of available panels from PanelApp.

        Args:
            base_url: PanelApp base URL

        Returns:
            List of panel information dictionaries
        """
        panels = []
        url: str | None = f"{base_url}/"

        while url:
            response = self._make_request(url)
            if not response:
                break

            panels.extend(response.get("results", []))
            url = response.get("next")  # Handle pagination

        logger.info(f"Found {len(panels)} panels from {base_url}")
        return panels

    def get_panel_genes(self, base_url: str, panel_id: int) -> list[dict[str, Any]]:
        """
        Get genes for a specific panel.

        Args:
            base_url: PanelApp base URL
            panel_id: Panel ID

        Returns:
            List of gene information dictionaries
        """
        genes = []
        url: str | None = f"{base_url}/{panel_id}/genes/"

        while url:
            response = self._make_request(url)
            if not response:
                break

            genes.extend(response.get("results", []))
            url = response.get("next")  # Handle pagination

        logger.info(f"Found {len(genes)} genes for panel {panel_id}")
        return genes

    def get_panel_info(self, base_url: str, panel_id: int) -> dict[str, Any] | None:
        """
        Get information about a specific panel.

        Args:
            base_url: PanelApp base URL
            panel_id: Panel ID

        Returns:
            Panel information dictionary or None
        """
        url = f"{base_url}/{panel_id}/"
        return self._make_request(url)


def fetch_panelapp_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch gene panel data from PanelApp APIs with caching support.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with PanelApp data
    """
    config_manager = ConfigManager(config)
    panelapp_config = config_manager.get_source_config("PanelApp")

    if not panelapp_config.get("enabled", True):
        logger.info("PanelApp data source is disabled")
        return pd.DataFrame()

    # Initialize cache manager
    cache_config = config.get("performance", {})
    cache_manager = CacheManager(
        cache_dir=cache_config.get("cache_dir", ".cache"),
        cache_ttl=cache_config.get("cache_ttl", 2592000),  # 30 days
        enabled=cache_config.get("enable_caching", True),
    )

    client = PanelAppClient(cache_manager=cache_manager)
    all_dataframes = []

    # Evidence level to score mapping
    evidence_scores = panelapp_config.get(
        "evidence_scores", {"Green": 1.0, "Amber": 0.5, "Red": 0.1}
    )

    # Process each PanelApp instance
    for panelapp_instance in panelapp_config.get("panels", []):
        instance_name = panelapp_instance.get("name", "PanelApp")
        base_url = panelapp_instance.get("base_url")

        if not base_url:
            logger.warning(f"No base URL specified for {instance_name}")
            continue

        logger.info(f"Fetching data from {instance_name}")

        # Get specified panels or all panels
        panels_to_fetch = panelapp_instance.get("panels", [])

        if not panels_to_fetch:
            # If no specific panels specified, get all panels
            logger.info(
                f"No specific panels configured, fetching all panels from {instance_name}"
            )
            all_panels = client.get_panel_list(base_url)
            panels_to_fetch = [
                {"id": panel["id"], "name": panel.get("name", f"Panel_{panel['id']}")}
                for panel in all_panels
            ]

        # Process each panel
        for panel_config in panels_to_fetch:
            panel_id = panel_config.get("id")
            panel_name = panel_config.get("name", f"Panel_{panel_id}")

            if not panel_id:
                logger.warning(f"No panel ID specified for {panel_name}")
                continue

            logger.info(f"Processing panel: {panel_name} (ID: {panel_id})")

            # Get panel genes
            panel_genes = client.get_panel_genes(base_url, panel_id)

            if not panel_genes:
                logger.warning(f"No genes found for panel {panel_name}")
                continue

            # Process genes
            genes = []
            scores = []
            details = []
            reported_names = []

            for gene_data in panel_genes:
                # Extract gene symbol
                gene_symbol = gene_data.get("gene_data", {}).get("gene_symbol", "")
                if not gene_symbol:
                    continue

                # Extract confidence level
                confidence = gene_data.get("confidence_level", "").strip()
                evidence_score = evidence_scores.get(confidence, 0.0)

                # Skip genes with zero evidence score
                if evidence_score == 0.0:
                    continue

                # Extract additional details
                mode_of_inheritance = gene_data.get("mode_of_inheritance", "")
                phenotypes = gene_data.get("phenotypes", [])
                publications = gene_data.get("publications", [])

                # Handle different phenotype formats
                if phenotypes and isinstance(phenotypes[0], dict):
                    # Old format: list of dictionaries
                    phenotype_names = [
                        p.get("phenotype", "") for p in phenotypes if p.get("phenotype")
                    ]
                else:
                    # New format: list of strings
                    phenotype_names = [str(p) for p in phenotypes if p]

                # Handle different publication formats
                if publications and isinstance(publications[0], dict):
                    # Old format: list of dictionaries
                    publication_ids = [
                        str(p.get("pmid", "")) for p in publications if p.get("pmid")
                    ]
                else:
                    # New format: list of strings
                    publication_ids = [str(p) for p in publications if p]

                detail_parts = [
                    f"Confidence:{confidence}",
                    f"MOI:{mode_of_inheritance}" if mode_of_inheritance else "",
                    (
                        f"Phenotypes:{';'.join(phenotype_names)}"
                        if phenotype_names
                        else ""
                    ),
                    f"PMIDs:{';'.join(publication_ids)}" if publication_ids else "",
                ]

                detail_string = "|".join([part for part in detail_parts if part])

                genes.append(gene_symbol)
                scores.append(evidence_score)
                details.append(detail_string)
                reported_names.append(gene_symbol)  # PanelApp uses standard symbols

            if genes:
                source_name = f"{instance_name}:{panel_name}"
                df = create_standard_dataframe(
                    genes=genes,
                    source_name=source_name,
                    evidence_scores=scores,
                    source_details=details,
                    gene_names_reported=reported_names,
                )
                all_dataframes.append(df)
                logger.info(f"Processed {len(genes)} genes from {source_name}")

    # Combine all dataframes
    if all_dataframes:
        combined_df = pd.concat(all_dataframes, ignore_index=True)
        logger.info(f"Fetched {len(combined_df)} total gene records from PanelApp")
        return combined_df
    else:
        logger.warning("No PanelApp data was successfully fetched")
        return pd.DataFrame()


def search_panelapp_panels(
    config: dict[str, Any], search_term: str
) -> list[dict[str, Any]]:
    """
    Search for panels in PanelApp by name or description.

    Args:
        config: Configuration dictionary
        search_term: Search term

    Returns:
        List of matching panels
    """
    config_manager = ConfigManager(config)
    panelapp_config = config_manager.get_source_config("PanelApp")

    if not panelapp_config.get("enabled", True):
        return []

    client = PanelAppClient()
    matching_panels = []

    # Search in each PanelApp instance
    for panelapp_instance in panelapp_config.get("panels", []):
        base_url = panelapp_instance.get("base_url")
        if not base_url:
            continue

        # Get all panels and filter by search term
        all_panels = client.get_panel_list(base_url)

        for panel in all_panels:
            panel_name = panel.get("name", "").lower()
            panel_description = panel.get("description", "").lower()

            if (
                search_term.lower() in panel_name
                or search_term.lower() in panel_description
            ):
                matching_panels.append(
                    {
                        "id": panel.get("id"),
                        "name": panel.get("name"),
                        "description": panel.get("description"),
                        "version": panel.get("version"),
                        "source": panelapp_instance.get("name", "PanelApp"),
                        "base_url": base_url,
                    }
                )

    return matching_panels
