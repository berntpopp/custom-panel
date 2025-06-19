"""
Ensembl client for gene annotation and genomic coordinate retrieval.

This module provides a client for interacting with the Ensembl REST API to
retrieve gene coordinates, transcript information, and other genomic data.
"""

import functools
import logging
import time
from typing import Any, Optional

import requests

logger = logging.getLogger(__name__)


class EnsemblClient:
    """Client for interacting with the Ensembl REST API."""

    BASE_URL = "https://rest.ensembl.org"

    def __init__(
        self, timeout: int = 30, max_retries: int = 3, retry_delay: float = 1.0
    ):
        """
        Initialize the Ensembl client.

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
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": "custom-panel/0.1.0",
            }
        )

    def _make_request(
        self, endpoint: str, method: str = "GET", data: dict[str, Any] | None = None
    ) -> dict[str, Any] | list[dict[str, Any]]:
        """
        Make a request to the Ensembl API with retry logic.

        Args:
            endpoint: API endpoint
            method: HTTP method (GET or POST)
            data: Request data for POST requests

        Returns:
            JSON response

        Raises:
            requests.RequestException: If request fails after retries
        """
        url = f"{self.BASE_URL}/{endpoint}"

        for attempt in range(self.max_retries + 1):
            try:
                if method.upper() == "POST":
                    response = self.session.post(url, json=data, timeout=self.timeout)
                else:
                    response = self.session.get(url, timeout=self.timeout)

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
    def get_gene_coordinates(
        self, gene_symbol: str, species: str = "human"
    ) -> Optional[dict[str, Any]]:
        """
        Get genomic coordinates for a gene symbol.

        Args:
            gene_symbol: Gene symbol
            species: Species name (default: "human")

        Returns:
            Dictionary with gene coordinates or None if not found
        """
        try:
            response = self._make_request(f"lookup/symbol/{species}/{gene_symbol}")
            if isinstance(response, dict):
                return {
                    "gene_id": response.get("id"),
                    "chromosome": response.get("seq_region_name"),
                    "start": response.get("start"),
                    "end": response.get("end"),
                    "strand": response.get("strand"),
                    "biotype": response.get("biotype"),
                    "description": response.get("description"),
                }
        except (requests.RequestException, ValueError):
            pass
        return None

    def get_genes_coordinates(
        self, gene_symbols: list[str], species: str = "human"
    ) -> dict[str, dict[str, Any] | None]:
        """
        Get genomic coordinates for multiple gene symbols using batch request.

        Args:
            gene_symbols: List of gene symbols
            species: Species name (default: "human")

        Returns:
            Dictionary mapping gene symbols to coordinate information
        """
        # Ensembl batch lookup
        try:
            data = {"symbols": gene_symbols}
            response = self._make_request(
                f"lookup/symbol/{species}", method="POST", data=data
            )

            result: dict[str, dict[str, Any] | None] = {}
            if isinstance(response, dict):
                for symbol in gene_symbols:
                    gene_data = response.get(symbol)
                    if gene_data:
                        result[symbol] = {
                            "gene_id": gene_data.get("id"),
                            "chromosome": gene_data.get("seq_region_name"),
                            "start": gene_data.get("start"),
                            "end": gene_data.get("end"),
                            "strand": gene_data.get("strand"),
                            "biotype": gene_data.get("biotype"),
                            "description": gene_data.get("description"),
                        }
                    else:
                        result[symbol] = None
            return result
        except requests.RequestException:
            # Fallback to individual requests
            return {
                symbol: self.get_gene_coordinates(symbol, species)
                for symbol in gene_symbols
            }

    @functools.lru_cache(maxsize=2000)  # noqa: B019
    def get_canonical_transcript(
        self, gene_id: str, species: str = "human"
    ) -> Optional[dict[str, Any]]:
        """
        Get the canonical transcript for a gene.

        Args:
            gene_id: Ensembl gene ID
            species: Species name (default: "human")

        Returns:
            Dictionary with canonical transcript information or None if not found
        """
        try:
            response = self._make_request(f"lookup/id/{gene_id}?expand=1")
            if isinstance(response, dict):
                transcripts = response.get("Transcript", [])
                # Find canonical transcript
                for transcript in transcripts:
                    if transcript.get("is_canonical") == 1:
                        return {
                            "transcript_id": transcript.get("id"),
                            "transcript_start": transcript.get("start"),
                            "transcript_end": transcript.get("end"),
                            "biotype": transcript.get("biotype"),
                            "length": transcript.get("length"),
                        }
                # If no canonical transcript found, return the first one
                if transcripts:
                    transcript = transcripts[0]
                    return {
                        "transcript_id": transcript.get("id"),
                        "transcript_start": transcript.get("start"),
                        "transcript_end": transcript.get("end"),
                        "biotype": transcript.get("biotype"),
                        "length": transcript.get("length"),
                    }
        except (requests.RequestException, ValueError):
            pass
        return None

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def get_mane_transcript(
        self, gene_symbol: str, species: str = "human"
    ) -> Optional[dict[str, Any]]:
        """
        Get MANE (Matched Annotation from NCBI and EBI) transcript information.

        Args:
            gene_symbol: Gene symbol
            species: Species name (default: "human")

        Returns:
            Dictionary with MANE transcript information or None if not found
        """
        # First get gene coordinates
        gene_info = self.get_gene_coordinates(gene_symbol, species)
        if not gene_info or not gene_info.get("gene_id"):
            return None

        # Get transcript information
        try:
            response = self._make_request(f"lookup/id/{gene_info['gene_id']}?expand=1")
            if isinstance(response, dict):
                transcripts = response.get("Transcript", [])
                for transcript in transcripts:
                    # Look for MANE select or MANE plus clinical transcripts
                    if any(
                        tag.get("value") in ["MANE Select", "MANE Plus Clinical"]
                        for tag in transcript.get("transcript_tags", [])
                    ):
                        return {
                            "transcript_id": transcript.get("id"),
                            "mane_type": next(
                                (
                                    tag.get("value")
                                    for tag in transcript.get("transcript_tags", [])
                                    if tag.get("value")
                                    in ["MANE Select", "MANE Plus Clinical"]
                                ),
                                None,
                            ),
                            "transcript_start": transcript.get("start"),
                            "transcript_end": transcript.get("end"),
                            "biotype": transcript.get("biotype"),
                            "length": transcript.get("length"),
                        }
        except (requests.RequestException, ValueError):
            pass
        return None

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def rsid_to_coordinates(
        self, rsid: str, species: str = "human"
    ) -> Optional[dict[str, Any]]:
        """
        Convert rsID to genomic coordinates.

        Args:
            rsid: dbSNP rsID (e.g., "rs1234567")
            species: Species name (default: "human")

        Returns:
            Dictionary with variant coordinates or None if not found
        """
        try:
            response = self._make_request(f"variation/{species}/{rsid}")
            if isinstance(response, dict):
                mappings = response.get("mappings", [])
                if mappings:
                    mapping = mappings[0]  # Take the first mapping
                    return {
                        "chromosome": mapping.get("seq_region_name"),
                        "start": mapping.get("start"),
                        "end": mapping.get("end"),
                        "strand": mapping.get("strand"),
                        "allele_string": mapping.get("allele_string"),
                        "assembly": mapping.get("coord_system"),
                    }
        except (requests.RequestException, ValueError):
            pass
        return None

    def calculate_gene_size(
        self, gene_symbol: str, species: str = "human"
    ) -> int | None:
        """
        Calculate the genomic size of a gene (end - start + 1).

        Args:
            gene_symbol: Gene symbol
            species: Species name (default: "human")

        Returns:
            Gene size in base pairs or None if not found
        """
        coords = self.get_gene_coordinates(gene_symbol, species)
        if coords and coords.get("start") and coords.get("end"):
            return coords["end"] - coords["start"] + 1
        return None

    def get_gene_annotation(
        self, gene_symbol: str, species: str = "human"
    ) -> dict[str, Any]:
        """
        Get comprehensive gene annotation information.

        Args:
            gene_symbol: Gene symbol
            species: Species name (default: "human")

        Returns:
            Dictionary with comprehensive gene annotation
        """
        result: dict[str, Any] = {
            "gene_symbol": gene_symbol,
            "coordinates": self.get_gene_coordinates(gene_symbol, species),
            "gene_size": None,
            "canonical_transcript": None,
            "mane_transcript": None,
        }

        # Calculate gene size
        if result["coordinates"]:
            result["gene_size"] = self.calculate_gene_size(gene_symbol, species)

            # Get canonical transcript
            coords = result["coordinates"]
            if coords and isinstance(coords, dict):
                gene_id = coords.get("gene_id")
                if gene_id:
                    result["canonical_transcript"] = self.get_canonical_transcript(
                        gene_id, species
                    )

        # Get MANE transcript
        result["mane_transcript"] = self.get_mane_transcript(gene_symbol, species)

        return result

    def get_cache_info(self) -> dict[str, Any]:
        """
        Get cache statistics for monitoring performance.

        Returns:
            Dictionary with cache statistics
        """
        return {
            "get_gene_coordinates": self.get_gene_coordinates.cache_info()._asdict(),
            "get_canonical_transcript": self.get_canonical_transcript.cache_info()._asdict(),
            "get_mane_transcript": self.get_mane_transcript.cache_info()._asdict(),
            "rsid_to_coordinates": self.rsid_to_coordinates.cache_info()._asdict(),
        }

    def clear_cache(self) -> None:
        """Clear all cached results."""
        self.get_gene_coordinates.cache_clear()
        self.get_canonical_transcript.cache_clear()
        self.get_mane_transcript.cache_clear()
        self.rsid_to_coordinates.cache_clear()
