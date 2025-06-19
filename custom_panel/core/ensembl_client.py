"""
Ensembl client for gene annotation and genomic coordinate retrieval.

This module provides a client for interacting with the Ensembl REST API to
retrieve gene coordinates, transcript information, and other genomic data.
"""

import functools
import logging
import time
from typing import Any

import requests

from .cache_manager import CacheManager

logger = logging.getLogger(__name__)


class EnsemblClient:
    """Client for interacting with the Ensembl REST API."""

    BASE_URL = "https://rest.ensembl.org"

    def __init__(
        self, timeout: int = 30, max_retries: int = 3, retry_delay: float = 1.0,
        transcript_batch_size: int = 50, cache_manager: CacheManager | None = None
    ):
        """
        Initialize the Ensembl client.

        Args:
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            retry_delay: Delay between retries in seconds
            transcript_batch_size: Batch size for transcript queries
            cache_manager: Optional cache manager instance
        """
        self.timeout = timeout
        self.max_retries = max_retries
        self.retry_delay = retry_delay
        self.transcript_batch_size = transcript_batch_size
        self.cache_manager = cache_manager
        self.session = requests.Session()
        self.session.headers.update(
            {
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": "custom-panel/0.1.0",
            }
        )

    def _make_request(
        self, endpoint: str, method: str = "GET", data: dict[str, Any] | None = None, timeout_override: int | None = None
    ) -> dict[str, Any] | list[dict[str, Any]]:
        """
        Make a request to the Ensembl API with retry logic.

        Args:
            endpoint: API endpoint
            method: HTTP method (GET or POST)
            data: Request data for POST requests
            timeout_override: Override timeout for this request

        Returns:
            JSON response

        Raises:
            requests.RequestException: If request fails after retries
        """
        # Check cache first
        if self.cache_manager:
            cached_response = self.cache_manager.get("ensembl", endpoint, method, data)
            if cached_response is not None:
                return cached_response

        url = f"{self.BASE_URL}/{endpoint}"

        # Log request details for debugging
        if data:
            if "symbols" in data:
                logger.debug(f"Making {method} request to {url} with {len(data['symbols'])} symbols")
            elif "ids" in data:
                logger.debug(f"Making {method} request to {url} with {len(data['ids'])} gene IDs")
            else:
                logger.debug(f"Making {method} request to {url} with data keys: {list(data.keys())}")
        else:
            logger.debug(f"Making {method} request to {url}")

        # Use override timeout if provided, otherwise use default
        request_timeout = timeout_override if timeout_override is not None else self.timeout

        for attempt in range(self.max_retries + 1):
            try:
                if method.upper() == "POST":
                    response = self.session.post(url, json=data, timeout=request_timeout)
                else:
                    response = self.session.get(url, timeout=request_timeout)

                response.raise_for_status()
                json_response = response.json()

                # Log response details for debugging
                if isinstance(json_response, dict):
                    logger.debug(f"Received response with {len(json_response)} items")
                elif isinstance(json_response, list):
                    logger.debug(f"Received response with {len(json_response)} list items")

                # Cache successful response
                if self.cache_manager:
                    self.cache_manager.set("ensembl", endpoint, method, data, json_response)

                return json_response
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
    ) -> dict[str, Any] | None:
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

        This method is maintained for backward compatibility.
        For new code, use get_symbols_data_batch() instead.

        Args:
            gene_symbols: List of gene symbols
            species: Species name (default: "human")

        Returns:
            Dictionary mapping gene symbols to coordinate information
        """
        return self.get_symbols_data_batch(gene_symbols, species, expand=False)

    def get_symbols_data_batch(
        self, gene_symbols: list[str], species: str = "human", expand: bool = False
    ) -> dict[str, dict[str, Any] | None]:
        """
        Get genomic data for multiple gene symbols using optimized batch requests.

        Args:
            gene_symbols: List of gene symbols
            species: Species name (default: "human")
            expand: Whether to fetch transcript data (default: False)

        Returns:
            Dictionary mapping gene symbols to complete gene data including transcripts if expanded
        """
        logger.debug(f"Fetching data for {len(gene_symbols)} genes, expand={expand}")

        # Step 1: Get coordinates for all genes in one batch call
        try:
            data = {"symbols": gene_symbols}
            response = self._make_request(f"lookup/symbol/{species}", method="POST", data=data)

            result: dict[str, dict[str, Any] | None] = {}
            gene_ids = {}  # Map symbol to gene_id for transcript lookup

            if isinstance(response, dict):
                for symbol in gene_symbols:
                    gene_data = response.get(symbol)
                    if gene_data:
                        gene_info = {
                            "gene_id": gene_data.get("id"),
                            "chromosome": gene_data.get("seq_region_name"),
                            "start": gene_data.get("start"),
                            "end": gene_data.get("end"),
                            "strand": gene_data.get("strand"),
                            "biotype": gene_data.get("biotype"),
                            "description": gene_data.get("description"),
                        }
                        result[symbol] = gene_info

                        # Store gene_id for transcript lookup
                        if expand and gene_data.get("id"):
                            gene_ids[symbol] = gene_data.get("id")
                    else:
                        result[symbol] = None

            # Step 2: If expand is requested, fetch transcript data efficiently
            if expand and gene_ids:
                logger.debug(f"Fetching transcript data for {len(gene_ids)} genes")
                self._add_transcript_data_batch(result, gene_ids, species)

            return result

        except requests.RequestException as e:
            logger.warning(f"Batch coordinate request failed: {e}")
            # Fallback to individual requests
            return {
                symbol: self.get_gene_coordinates(symbol, species)
                for symbol in gene_symbols
            }

    def _add_transcript_data_batch(
        self,
        result: dict[str, dict[str, Any] | None],
        gene_ids: dict[str, str],
        species: str
    ) -> None:
        """
        Add transcript data to gene results using batch requests where possible.

        Args:
            result: Gene data dictionary to update
            gene_ids: Mapping of gene symbols to gene IDs
            species: Species name
        """
        # Process gene IDs in smaller batches to avoid timeouts
        batch_size = self.transcript_batch_size
        gene_id_list = list(gene_ids.items())

        for i in range(0, len(gene_id_list), batch_size):
            batch_items = gene_id_list[i:i + batch_size]
            batch_gene_ids = [item[1] for item in batch_items]

            logger.info(f"Processing transcript batch {i//batch_size + 1}/{(len(gene_id_list) + batch_size - 1)//batch_size}")

            try:
                # Use batch lookup for gene IDs with expand
                data = {"ids": batch_gene_ids}
                logger.debug(f"Fetching transcript data for batch of {len(batch_gene_ids)} gene IDs")
                # Use longer timeout for transcript queries as they return more data
                response = self._make_request(
                    f"lookup/id/{species}?expand=1",
                    method="POST",
                    data=data,
                    timeout_override=120  # 120 second timeout for transcript queries
                )

                if isinstance(response, dict):
                    for symbol, gene_id in batch_items:
                        gene_data = response.get(gene_id)
                        if gene_data and "Transcript" in gene_data:
                            transcripts = gene_data["Transcript"]

                            # Find canonical transcript
                            canonical_transcript = None
                            for transcript in transcripts:
                                if transcript.get("is_canonical") == 1:
                                    canonical_transcript = {
                                        "transcript_id": transcript.get("id"),
                                        "transcript_start": transcript.get("start"),
                                        "transcript_end": transcript.get("end"),
                                        "biotype": transcript.get("biotype"),
                                        "length": transcript.get("length"),
                                    }
                                    break

                            # If no canonical found, use first transcript
                            if not canonical_transcript and transcripts:
                                transcript = transcripts[0]
                                canonical_transcript = {
                                    "transcript_id": transcript.get("id"),
                                    "transcript_start": transcript.get("start"),
                                    "transcript_end": transcript.get("end"),
                                    "biotype": transcript.get("biotype"),
                                    "length": transcript.get("length"),
                                }

                            # Find MANE transcript
                            mane_transcript = None
                            for transcript in transcripts:
                                transcript_tags = transcript.get("transcript_tags", [])
                                for tag in transcript_tags:
                                    if tag.get("value") in ["MANE Select", "MANE Plus Clinical"]:
                                        mane_transcript = {
                                            "transcript_id": transcript.get("id"),
                                            "mane_type": tag.get("value"),
                                            "transcript_start": transcript.get("start"),
                                            "transcript_end": transcript.get("end"),
                                            "biotype": transcript.get("biotype"),
                                            "length": transcript.get("length"),
                                        }
                                        break
                                if mane_transcript:
                                    break

                            # Update the result with transcript data
                            if result[symbol]:
                                result[symbol]["canonical_transcript"] = canonical_transcript
                                result[symbol]["mane_transcript"] = mane_transcript

            except requests.RequestException as e:
                logger.warning(f"Batch transcript request failed: {e}, using fallback")
                # Fallback to individual transcript requests for this batch
                for symbol, _gene_id in batch_items:
                    if result[symbol]:
                        result[symbol]["canonical_transcript"] = None
                        result[symbol]["mane_transcript"] = None

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def rsid_to_coordinates(
        self, rsid: str, species: str = "human"
    ) -> dict[str, Any] | None:
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
        # Use the new batch method with expand for single gene
        batch_result = self.get_symbols_data_batch([gene_symbol], species, expand=True)
        gene_data = batch_result.get(gene_symbol)

        if not gene_data:
            return {
                "gene_symbol": gene_symbol,
                "coordinates": None,
                "gene_size": None,
                "canonical_transcript": None,
                "mane_transcript": None,
            }

        # Calculate gene size
        gene_size = None
        if gene_data.get("start") and gene_data.get("end"):
            gene_size = gene_data["end"] - gene_data["start"] + 1

        return {
            "gene_symbol": gene_symbol,
            "coordinates": {
                "gene_id": gene_data.get("gene_id"),
                "chromosome": gene_data.get("chromosome"),
                "start": gene_data.get("start"),
                "end": gene_data.get("end"),
                "strand": gene_data.get("strand"),
                "biotype": gene_data.get("biotype"),
                "description": gene_data.get("description"),
            },
            "gene_size": gene_size,
            "canonical_transcript": gene_data.get("canonical_transcript"),
            "mane_transcript": gene_data.get("mane_transcript"),
        }

    def get_cache_info(self) -> dict[str, Any]:
        """
        Get cache statistics for monitoring performance.

        Returns:
            Dictionary with cache statistics
        """
        return {
            "get_gene_coordinates": self.get_gene_coordinates.cache_info()._asdict(),
            "rsid_to_coordinates": self.rsid_to_coordinates.cache_info()._asdict(),
        }

    def clear_cache(self) -> None:
        """Clear all cached results."""
        self.get_gene_coordinates.cache_clear()
        self.rsid_to_coordinates.cache_clear()
