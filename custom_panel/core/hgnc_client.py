"""
HGNC (HUGO Gene Nomenclature Committee) client for gene symbol standardization.

This module provides a client for interacting with the HGNC REST API to
standardize gene symbols and retrieve gene information.
"""

import functools
import logging
import time
from typing import Any

import requests

logger = logging.getLogger(__name__)


class HGNCClient:
    """Client for interacting with the HGNC REST API."""

    BASE_URL = "https://rest.genenames.org"

    def __init__(
        self, timeout: int = 30, max_retries: int = 3, retry_delay: float = 1.0
    ):
        """
        Initialize the HGNC client.

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
            {"Accept": "application/json", "User-Agent": "custom-panel/0.1.0"}
        )

    def _make_request(
        self,
        endpoint: str,
        params: dict[str, Any] | None = None,
        method: str = "GET",
        data: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """
        Make a request to the HGNC API with retry logic.

        Args:
            endpoint: API endpoint
            params: Query parameters
            method: HTTP method (GET or POST)
            data: Request data for POST requests

        Returns:
            JSON response as dictionary

        Raises:
            requests.RequestException: If request fails after retries
        """
        url = f"{self.BASE_URL}/{endpoint}"

        for attempt in range(self.max_retries + 1):
            try:
                if method.upper() == "POST":
                    response = self.session.post(
                        url, json=data, params=params, timeout=self.timeout
                    )
                else:
                    response = self.session.get(
                        url, params=params, timeout=self.timeout
                    )
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

    @functools.lru_cache(maxsize=10000)  # noqa: B019
    def symbol_to_hgnc_id(self, symbol: str) -> str | None:
        """
        Convert a gene symbol to HGNC ID.

        Args:
            symbol: Gene symbol

        Returns:
            HGNC ID (e.g., "HGNC:5") or None if not found
        """
        try:
            response = self._make_request("search/symbol", {"symbol": symbol})
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0].get("hgnc_id")
        except (requests.RequestException, ValueError):
            pass
        return None

    @functools.lru_cache(maxsize=10000)  # noqa: B019
    def hgnc_id_to_symbol(self, hgnc_id: str) -> str | None:
        """
        Convert HGNC ID to approved gene symbol.

        Args:
            hgnc_id: HGNC ID (e.g., "HGNC:5")

        Returns:
            Approved gene symbol or None if not found
        """
        try:
            response = self._make_request("search/hgnc_id", {"hgnc_id": hgnc_id})
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0].get("symbol")
        except (requests.RequestException, ValueError):
            pass
        return None

    @functools.lru_cache(maxsize=10000)  # noqa: B019
    def get_gene_info(self, symbol: str) -> dict[str, Any] | None:
        """
        Get comprehensive gene information for a symbol.

        Args:
            symbol: Gene symbol

        Returns:
            Dictionary with gene information or None if not found
        """
        try:
            response = self._make_request("search/symbol", {"symbol": symbol})
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0]
        except (requests.RequestException, ValueError):
            pass
        return None

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def standardize_symbol(self, symbol: str) -> str:
        """
        Standardize a gene symbol to the approved HGNC symbol.

        This method tries multiple approaches:
        1. Direct symbol lookup
        2. Previous symbol lookup
        3. Alias symbol lookup

        Args:
            symbol: Input gene symbol

        Returns:
            Standardized HGNC symbol or original symbol if not found
        """
        symbol = symbol.strip().upper()

        # Try direct symbol lookup
        gene_info = self.get_gene_info(symbol)
        if gene_info:
            return gene_info.get("symbol", symbol)

        # Try previous symbol lookup
        try:
            response = self._make_request("search/prev_symbol", {"prev_symbol": symbol})
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0].get("symbol", symbol)
        except (requests.RequestException, ValueError):
            pass

        # Try alias symbol lookup
        try:
            response = self._make_request(
                "search/alias_symbol", {"alias_symbol": symbol}
            )
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0].get("symbol", symbol)
        except (requests.RequestException, ValueError):
            pass

        logger.warning(f"Could not standardize gene symbol: {symbol}")
        return symbol

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def standardize_symbols_batch(self, symbols: tuple[str, ...]) -> dict[str, str]:
        """
        Standardize multiple gene symbols using HGNC batch API.

        Args:
            symbols: Tuple of gene symbols (tuple for caching compatibility)

        Returns:
            Dictionary mapping original symbols to standardized symbols
        """
        if not symbols:
            return {}

        # Keep track of original input symbols and their normalized versions
        original_symbols = list(symbols)
        normalized_symbols = [symbol.strip().upper() for symbol in symbols]
        symbol_map = dict(zip(normalized_symbols, original_symbols, strict=True))

        # Try batch direct symbol lookup
        result = {}
        try:
            # Use the HGNC batch search endpoint
            batch_data = {"symbols": normalized_symbols}
            response = self._make_request("search", method="POST", data=batch_data)

            # Parse response - HGNC batch endpoint returns a list of docs
            docs = response.get("response", {}).get("docs", [])

            # Create mapping from found symbols
            found_symbols = set()
            for doc in docs:
                standardized_symbol = doc.get("symbol")
                if standardized_symbol:
                    # Find which input symbol this matches (case-insensitive)
                    for normalized_input in normalized_symbols:
                        if normalized_input == standardized_symbol.upper():
                            original_input = symbol_map[normalized_input]
                            result[original_input] = standardized_symbol
                            found_symbols.add(normalized_input)
                            break

                    # Also check aliases and previous symbols
                    alias_symbols = doc.get("alias_symbol", [])
                    if isinstance(alias_symbols, list):
                        for alias in alias_symbols:
                            for normalized_input in normalized_symbols:
                                if normalized_input == alias.upper():
                                    original_input = symbol_map[normalized_input]
                                    result[original_input] = standardized_symbol
                                    found_symbols.add(normalized_input)
                                    break

                    prev_symbols = doc.get("prev_symbol", [])
                    if isinstance(prev_symbols, list):
                        for prev in prev_symbols:
                            for normalized_input in normalized_symbols:
                                if normalized_input == prev.upper():
                                    original_input = symbol_map[normalized_input]
                                    result[original_input] = standardized_symbol
                                    found_symbols.add(normalized_input)
                                    break

        except requests.RequestException as e:
            logger.warning(
                f"Batch symbol standardization failed: {e}, falling back to individual requests"
            )
            # Fall back to individual requests
            for original_symbol in original_symbols:
                result[original_symbol] = self.standardize_symbol(original_symbol)
            return result

        # For symbols not found in batch response, use original symbol
        for normalized_input in normalized_symbols:
            if normalized_input not in found_symbols:
                original_input = symbol_map[normalized_input]
                result[original_input] = original_input

        logger.info(
            f"Batch standardized {len([k for k, v in result.items() if k != v])} out of {len(symbols)} symbols"
        )
        return result

    def standardize_symbols(self, symbols: list[str]) -> dict[str, str]:
        """
        Standardize multiple gene symbols.

        Args:
            symbols: List of gene symbols

        Returns:
            Dictionary mapping original symbols to standardized symbols
        """
        # Convert to tuple for caching compatibility and call batch method
        return self.standardize_symbols_batch(tuple(symbols))

    def get_cache_info(self) -> dict[str, Any]:
        """
        Get cache statistics for monitoring performance.

        Returns:
            Dictionary with cache statistics
        """
        return {
            "symbol_to_hgnc_id": self.symbol_to_hgnc_id.cache_info()._asdict(),
            "hgnc_id_to_symbol": self.hgnc_id_to_symbol.cache_info()._asdict(),
            "get_gene_info": self.get_gene_info.cache_info()._asdict(),
            "standardize_symbol": self.standardize_symbol.cache_info()._asdict(),
            "standardize_symbols_batch": self.standardize_symbols_batch.cache_info()._asdict(),
        }

    def clear_cache(self) -> None:
        """Clear all cached results."""
        self.symbol_to_hgnc_id.cache_clear()
        self.hgnc_id_to_symbol.cache_clear()
        self.get_gene_info.cache_clear()
        self.standardize_symbol.cache_clear()
        self.standardize_symbols_batch.cache_clear()
