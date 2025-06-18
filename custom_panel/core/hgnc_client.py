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
        try:
            response = self._make_request(f"search/symbol/{symbol}", method="GET")
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0].get("symbol", symbol)
        except (requests.RequestException, ValueError):
            pass

        # Try previous symbol lookup
        try:
            response = self._make_request(f"search/prev_symbol/{symbol}", method="GET")
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0].get("symbol", symbol)
        except (requests.RequestException, ValueError):
            pass

        # Try alias symbol lookup
        try:
            response = self._make_request(f"search/alias_symbol/{symbol}", method="GET")
            docs = response.get("response", {}).get("docs", [])
            if docs:
                return docs[0].get("symbol", symbol)
        except (requests.RequestException, ValueError):
            pass

        logger.warning(f"Could not standardize gene symbol: {symbol}")
        return symbol

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def standardize_symbols_batch(
        self, symbols: tuple[str, ...]
    ) -> dict[str, dict[str, str | None]]:
        """
        Standardize multiple gene symbols using HGNC batch API.

        This method now constructs a proper GET request with an 'OR' query.

        Args:
            symbols: Tuple of gene symbols (tuple for caching compatibility)

        Returns:
            Dictionary mapping original symbols to dict containing approved_symbol and hgnc_id
        """
        if not symbols:
            return {}

        original_symbols = list(symbols)
        # The result dictionary will map original input symbols to their standardized info.
        result = {
            symbol: {"approved_symbol": symbol, "hgnc_id": None}
            for symbol in original_symbols
        }

        try:
            # Construct the query for the correct endpoint: "GENE1+OR+GENE2+OR+..."
            query_str = "+OR+".join([s.upper() for s in original_symbols])
            endpoint = f"search/symbol/{query_str}"

            # Use the correct search/symbol endpoint with fields parameter
            params = {"fields": "symbol,hgnc_id"}
            response = self._make_request(endpoint, params=params, method="GET")

            docs = response.get("response", {}).get("docs", [])

            # Just match by symbol name directly
            found_symbols = set()
            for doc in docs:
                approved_symbol = doc.get("symbol")
                hgnc_id = doc.get("hgnc_id")
                if not approved_symbol:
                    continue

                # Check if this approved symbol matches any input symbol (case-insensitive)
                for original_symbol in original_symbols:
                    if original_symbol.upper() == approved_symbol.upper():
                        result[original_symbol] = {
                            "approved_symbol": approved_symbol,
                            "hgnc_id": hgnc_id,
                        }
                        found_symbols.add(original_symbol)
                        break

            # Log batch results
            exact_matches = len(found_symbols)
            need_postprocess = len(original_symbols) - exact_matches
            logger.debug(
                f"HGNC batch API results: {len(original_symbols)} symbols submitted, {exact_matches} exact matches found, {need_postprocess} need further processing"
            )

            if need_postprocess > 0 and logger.isEnabledFor(logging.DEBUG):
                unmatched = [s for s in original_symbols if s not in found_symbols]
                logger.debug(
                    f"Symbols requiring alias/prev_symbol search: {unmatched[:10]}{'...' if len(unmatched) > 10 else ''}"
                )

            # For symbols not found in batch, try individual lookups with aliases/prev symbols
            postprocess_fixed = 0
            for original_symbol in original_symbols:
                if original_symbol not in found_symbols:
                    # Try individual standardization (includes alias and prev symbol lookups)
                    standardized = self.standardize_symbol(original_symbol)
                    if standardized.upper() != original_symbol.upper():
                        # Symbol was standardized, try to get its HGNC ID
                        gene_info = self.get_gene_info(standardized)
                        hgnc_id = gene_info.get("hgnc_id") if gene_info else None
                        result[original_symbol] = {
                            "approved_symbol": standardized,
                            "hgnc_id": hgnc_id,
                        }
                        postprocess_fixed += 1
                    else:
                        # Symbol not changed, keep the default
                        result[original_symbol] = {
                            "approved_symbol": standardized,
                            "hgnc_id": None,
                        }

            if need_postprocess > 0:
                logger.debug(
                    f"Alias/prev_symbol search results: {postprocess_fixed} out of {need_postprocess} symbols were successfully resolved"
                )

        except requests.RequestException as e:
            logger.warning(
                f"Batch symbol standardization failed: {e}. Falling back to individual lookups."
            )
            # Fallback to individual, slower lookups if batch fails
            for original_symbol in original_symbols:
                standardized = self.standardize_symbol(original_symbol)
                gene_info = self.get_gene_info(standardized)
                hgnc_id = gene_info.get("hgnc_id") if gene_info else None
                result[original_symbol] = {
                    "approved_symbol": standardized,
                    "hgnc_id": hgnc_id,
                }

        changed_count = sum(
            1
            for k, v in result.items()
            if v["approved_symbol"] is not None
            and k.upper() != v["approved_symbol"].upper()
        )
        logger.debug(
            f"HGNC batch complete: {len(symbols)} symbols processed, {changed_count} standardized to different symbols"
        )

        if changed_count > 0 and logger.isEnabledFor(logging.DEBUG):
            changes = [
                (orig, info["approved_symbol"])
                for orig, info in result.items()
                if info["approved_symbol"] is not None
                and orig.upper() != info["approved_symbol"].upper()
            ][:5]
            logger.debug(f"Example standardizations: {changes}")

        return result

    def standardize_symbols(
        self, symbols: list[str]
    ) -> dict[str, dict[str, str | None]]:
        """
        Standardize multiple gene symbols.

        Args:
            symbols: List of gene symbols

        Returns:
            Dictionary mapping original symbols to dict containing approved_symbol and hgnc_id
        """
        # For large lists, split into smaller batches to avoid URL length limits
        batch_size = 100  # Batch size for URL length limits
        result = {}

        for i in range(0, len(symbols), batch_size):
            batch = symbols[i : i + batch_size]
            batch_result = self.standardize_symbols_batch(tuple(batch))
            result.update(batch_result)

        return result

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
