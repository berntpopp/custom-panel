"""
HPO (Human Phenotype Ontology) client for phenotype and gene data retrieval.

This module provides a client for interacting with the HPO API to retrieve
phenotype terms, gene associations, and ontology relationships.
"""

import functools
import logging
import time
from typing import Any

import requests

logger = logging.getLogger(__name__)


class HPOClient:
    """Client for interacting with the HPO API."""

    BASE_URL = "'https://ontology.jax.org/api"

    def __init__(
        self, timeout: int = 30, max_retries: int = 3, retry_delay: float = 1.0
    ):
        """
        Initialize the HPO client.

        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            retry_delay: Delay between retries in seconds
        """
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.session = requests.Session()
        self.session.headers.update(
            {
                "Accept": "application/json",
                "User-Agent": "custom-panel/0.1.0",
            }
        )

    def _make_request(
        self, endpoint: str, params: dict[str, Any] | None = None
    ) -> dict[str, Any] | list[dict[str, Any]]:
        """
        Make a request to the HPO API with retry logic.

        Args:
            endpoint: API endpoint
            params: Query parameters

        Returns:
            JSON response

        Raises:
            requests.RequestException: If request fails after retries
        """
        url = f"{self.BASE_URL}/{endpoint}"

        for attempt in range(self.max_retries + 1):
            try:
                response = self.session.get(url, params=params, timeout=self.timeout)
                response.raise_for_status()
                return response.json()
            except (requests.RequestException, ValueError) as e:
                if attempt == self.max_retries:
                    logger.error(
                        f"Failed to fetch {url} after {self.max_retries} retries: {e}"
                    )
                    raise
                logger.warning(
                    f"Request failed (attempt {attempt + 1}/{self.max_retries + 1}): {e}"
                )
                time.sleep(self.retry_delay * (2**attempt))  # Exponential backoff

        # This should never be reached, but satisfies type checker
        raise requests.RequestException("Unexpected error in request loop")

    @functools.lru_cache(maxsize=5000)  # noqa: B019
    def get_term_info(self, hpo_id: str) -> dict[str, Any] | None:
        """
        Get information about an HPO term.

        Args:
            hpo_id: HPO term ID (e.g., "HP:0002664")

        Returns:
            Dictionary with term information or None if not found
        """
        try:
            response = self._make_request(f"hpo/term/{hpo_id}")
            if isinstance(response, dict):
                return {
                    "id": response.get("id"),
                    "name": response.get("name"),
                    "definition": response.get("definition"),
                    "synonyms": response.get("synonyms", []),
                    "parents": response.get("parents", []),
                    "children": response.get("children", []),
                    "genes": response.get("genes", []),
                }
        except (requests.RequestException, ValueError):
            pass
        return None

    @functools.lru_cache(maxsize=2000)  # noqa: B019
    def search_terms(self, query: str, limit: int = 100) -> list[dict[str, Any]]:
        """
        Search for HPO terms by name or synonym.

        Args:
            query: Search query
            limit: Maximum number of results

        Returns:
            List of matching HPO terms
        """
        try:
            params = {"q": query, "max": limit}
            response = self._make_request("hp/search", params=params)
            if isinstance(response, dict):
                terms = response.get("terms", [])
                return [
                    {
                        "id": term.get("id"),
                        "name": term.get("name"),
                        "definition": term.get("definition"),
                        "synonyms": term.get("synonyms", []),
                    }
                    for term in terms
                ]
        except (requests.RequestException, ValueError):
            pass
        return []

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def get_term_genes(self, hpo_id: str) -> list[dict[str, Any]]:
        """
        Get genes associated with an HPO term.

        Args:
            hpo_id: HPO term ID

        Returns:
            List of associated genes
        """
        try:
            response = self._make_request(f"hpo/term/{hpo_id}/genes")
            if isinstance(response, dict):
                genes = response.get("genes", [])
                return [
                    {
                        "gene_symbol": gene.get("gene_symbol"),
                        "entrez_id": gene.get("entrez_id"),
                        "disease_id": gene.get("disease_id"),
                        "disease_name": gene.get("disease_name"),
                        "frequency": gene.get("frequency"),
                        "onset": gene.get("onset"),
                        "evidence": gene.get("evidence"),
                    }
                    for gene in genes
                ]
        except (requests.RequestException, ValueError):
            pass
        return []

    def get_descendant_terms(
        self, hpo_id: str, max_depth: int = 10, include_self: bool = True
    ) -> set[str]:
        """
        Get all descendant terms of a given HPO term recursively.

        Args:
            hpo_id: HPO term ID
            max_depth: Maximum recursion depth
            include_self: Whether to include the original term

        Returns:
            Set of HPO term IDs including descendants
        """
        descendants = set()
        if include_self:
            descendants.add(hpo_id)

        def _collect_descendants(term_id: str, current_depth: int) -> None:
            if current_depth >= max_depth:
                return

            term_info = self.get_term_info(term_id)
            if term_info and "children" in term_info:
                children = term_info.get("children", [])
                for child in children:
                    child_id = child.get("id") if isinstance(child, dict) else child
                    if child_id and child_id not in descendants:
                        descendants.add(child_id)
                        _collect_descendants(child_id, current_depth + 1)

        _collect_descendants(hpo_id, 0)
        return descendants

    def get_genes_for_phenotype_hierarchy(
        self, root_hpo_id: str, max_depth: int = 10
    ) -> dict[str, dict[str, Any]]:
        """
        Get all genes associated with a phenotype and its descendants.

        Args:
            root_hpo_id: Root HPO term ID
            max_depth: Maximum hierarchy depth to traverse

        Returns:
            Dictionary mapping gene symbols to gene information
        """
        # Get all descendant terms
        all_terms = self.get_descendant_terms(
            root_hpo_id, max_depth=max_depth, include_self=True
        )

        logger.info(f"Found {len(all_terms)} terms in hierarchy for {root_hpo_id}")

        # Collect genes from all terms
        genes_dict = {}
        for term_id in all_terms:
            term_genes = self.get_term_genes(term_id)
            for gene in term_genes:
                gene_symbol = gene.get("gene_symbol")
                if gene_symbol:
                    if gene_symbol not in genes_dict:
                        genes_dict[gene_symbol] = {
                            "gene_symbol": gene_symbol,
                            "entrez_id": gene.get("entrez_id"),
                            "hpo_terms": [],
                            "diseases": set(),
                            "frequencies": [],
                            "evidence_sources": set(),
                        }

                    # Add HPO term association
                    genes_dict[gene_symbol]["hpo_terms"].append(
                        {
                            "hpo_id": term_id,
                            "disease_id": gene.get("disease_id"),
                            "disease_name": gene.get("disease_name"),
                            "frequency": gene.get("frequency"),
                            "onset": gene.get("onset"),
                            "evidence": gene.get("evidence"),
                        }
                    )

                    # Collect unique diseases and evidence
                    if gene.get("disease_name"):
                        genes_dict[gene_symbol]["diseases"].add(
                            gene.get("disease_name")
                        )
                    if gene.get("frequency"):
                        genes_dict[gene_symbol]["frequencies"].append(
                            gene.get("frequency")
                        )
                    if gene.get("evidence"):
                        genes_dict[gene_symbol]["evidence_sources"].add(
                            gene.get("evidence")
                        )

        # Convert sets to lists for JSON serialization
        for gene_data in genes_dict.values():
            gene_data["diseases"] = list(gene_data["diseases"])
            gene_data["evidence_sources"] = list(gene_data["evidence_sources"])

        logger.info(f"Found {len(genes_dict)} unique genes for hierarchy {root_hpo_id}")
        return genes_dict

    @functools.lru_cache(maxsize=100)  # noqa: B019
    def find_neoplasm_terms(self) -> list[dict[str, Any]]:
        """
        Find HPO terms related to neoplasms/cancer.

        Returns:
            List of neoplasm-related HPO terms
        """
        neoplasm_queries = [
            "neoplasm",
            "tumor",
            "tumour",
            "cancer",
            "carcinoma",
            "sarcoma",
            "lymphoma",
            "leukemia",
            "leukaemia",
            "melanoma",
            "glioma",
            "blastoma",
            "malignancy",
            "malignant",
        ]

        all_terms = []
        seen_ids = set()

        for query in neoplasm_queries:
            terms = self.search_terms(query, limit=50)
            for term in terms:
                term_id = term.get("id")
                if term_id and term_id not in seen_ids:
                    all_terms.append(term)
                    seen_ids.add(term_id)

        logger.info(f"Found {len(all_terms)} neoplasm-related HPO terms")
        return all_terms

    def get_omim_disease_genes(self, omim_id: str) -> list[dict[str, Any]]:
        """
        Get genes associated with an OMIM disease (if available through HPO).

        Args:
            omim_id: OMIM disease ID (e.g., "OMIM:114480")

        Returns:
            List of associated genes
        """
        try:
            # HPO API may have OMIM disease associations
            response = self._make_request(f"hpo/disease/{omim_id}")
            if isinstance(response, dict):
                genes = response.get("genes", [])
                return [
                    {
                        "gene_symbol": gene.get("gene_symbol"),
                        "entrez_id": gene.get("entrez_id"),
                        "hpo_terms": gene.get("hpo_terms", []),
                    }
                    for gene in genes
                ]
        except (requests.RequestException, ValueError):
            pass
        return []

    def get_cache_info(self) -> dict[str, Any]:
        """
        Get cache statistics for monitoring performance.

        Returns:
            Dictionary with cache statistics
        """
        return {
            "get_term_info": self.get_term_info.cache_info()._asdict(),
            "search_terms": self.search_terms.cache_info()._asdict(),
            "get_term_genes": self.get_term_genes.cache_info()._asdict(),
            "find_neoplasm_terms": self.find_neoplasm_terms.cache_info()._asdict(),
        }

    def clear_cache(self) -> None:
        """Clear all cached results."""
        self.get_term_info.cache_clear()
        self.search_terms.cache_clear()
        self.get_term_genes.cache_clear()
        self.find_neoplasm_terms.cache_clear()
