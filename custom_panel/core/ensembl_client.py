"""
Ensembl client for gene annotation and genomic coordinate retrieval.

This module provides a client for interacting with the Ensembl REST API to
retrieve gene coordinates, transcript information, and other genomic data.
"""

import functools
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import requests

from .cache_manager import CacheManager

logger = logging.getLogger(__name__)


class EnsemblClient:
    """Client for interacting with the Ensembl REST API."""

    BASE_URL = "https://rest.ensembl.org"

    def __init__(
        self,
        timeout: int = 30,
        max_retries: int = 3,
        retry_delay: float = 1.0,
        transcript_batch_size: int = 50,
        cache_manager: CacheManager | None = None,
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
            },
        )

    def _make_request(
        self,
        endpoint: str,
        method: str = "GET",
        data: dict[str, Any] | None = None,
        timeout_override: int | None = None,
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
                logger.debug(
                    f"Making {method} request to {url} with {len(data['symbols'])} symbols",
                )
            elif "ids" in data:
                logger.debug(
                    f"Making {method} request to {url} with {len(data['ids'])} gene IDs",
                )
            else:
                logger.debug(
                    f"Making {method} request to {url} with data keys: {list(data.keys())}",
                )
        else:
            logger.debug(f"Making {method} request to {url}")

        # Use override timeout if provided, otherwise use default
        request_timeout = (
            timeout_override if timeout_override is not None else self.timeout
        )

        for attempt in range(self.max_retries + 1):
            try:
                if method.upper() == "POST":
                    response = self.session.post(
                        url,
                        json=data,
                        timeout=request_timeout,
                    )
                else:
                    response = self.session.get(url, timeout=request_timeout)

                response.raise_for_status()
                json_response = response.json()

                # Log response details for debugging
                if isinstance(json_response, dict):
                    logger.debug(f"Received response with {len(json_response)} items")
                elif isinstance(json_response, list):
                    logger.debug(
                        f"Received response with {len(json_response)} list items",
                    )

                # Cache successful response
                if self.cache_manager:
                    self.cache_manager.set(
                        "ensembl",
                        endpoint,
                        method,
                        data,
                        json_response,
                    )

                return json_response
            except (requests.RequestException, ValueError) as e:
                if attempt == self.max_retries:
                    logger.error(
                        f"Failed to fetch {url} after {self.max_retries} retries: {e}",
                    )
                    raise
                logger.warning(
                    f"Request failed (attempt {attempt + 1}/{self.max_retries + 1}): {e}",
                )
                time.sleep(self.retry_delay * (2**attempt))  # Exponential backoff

        # This should never be reached, but satisfies type checker
        raise requests.RequestException("Unexpected error in request loop")

    @functools.lru_cache(maxsize=5000)  # noqa: B019
    def get_gene_coordinates(
        self,
        gene_symbol: str,
        species: str = "human",
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
        self,
        gene_symbols: list[str],
        species: str = "human",
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
        self,
        gene_symbols: list[str],
        species: str = "human",
        expand: bool = False,
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
            response = self._make_request(
                f"lookup/symbol/{species}",
                method="POST",
                data=data,
            )

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
        species: str,
    ) -> None:
        """
        Add transcript data to gene results using parallel batch requests.

        Args:
            result: Gene data dictionary to update
            gene_ids: Mapping of gene symbols to gene IDs
            species: Species name
        """
        # Process gene IDs in smaller batches to avoid timeouts
        batch_size = self.transcript_batch_size
        gene_id_list = list(gene_ids.items())

        # Split into batches
        batches = [
            gene_id_list[i : i + batch_size]
            for i in range(0, len(gene_id_list), batch_size)
        ]

        if len(batches) == 1:
            # Single batch - process directly
            self._process_transcript_batch(batches[0], result, species)
            return

        # Multiple batches - use parallel processing with limited workers
        # Use fewer workers for transcript requests as they're memory-intensive
        max_transcript_workers = min(3, len(batches))
        logger.info(
            f"Processing {len(batches)} transcript batches in parallel (workers={max_transcript_workers})",
        )

        with ThreadPoolExecutor(max_workers=max_transcript_workers) as executor:
            # Submit all transcript batch jobs
            future_to_batch = {
                executor.submit(
                    self._process_transcript_batch,
                    batch_items,
                    result,
                    species,
                ): i
                for i, batch_items in enumerate(batches)
            }

            # Collect results as they complete
            completed = 0
            for future in as_completed(future_to_batch):
                batch_num = future_to_batch[future]
                try:
                    future.result()  # This will raise if there was an exception
                    completed += 1
                    logger.info(
                        f"✓ Completed transcript batch {completed}/{len(batches)}",
                    )
                except Exception as e:
                    logger.error(f"✗ Transcript batch {batch_num + 1} failed: {e}")
                    # The batch method already handles fallbacks internally

    def _process_transcript_batch(
        self,
        batch_items: list[tuple[str, str]],
        result: dict[str, dict[str, Any] | None],
        species: str,
    ) -> None:
        """
        Process a single transcript batch.

        Args:
            batch_items: List of (symbol, gene_id) tuples
            result: Gene data dictionary to update
            species: Species name
        """
        batch_gene_ids = [item[1] for item in batch_items]

        try:
            # Use batch lookup for gene IDs with expand and MANE transcript data
            data = {"ids": batch_gene_ids}
            logger.debug(
                f"Fetching transcript data for batch of {len(batch_gene_ids)} gene IDs",
            )
            # Use longer timeout for transcript queries as they return more data
            response = self._make_request(
                f"lookup/id/{species}?expand=1&mane=1",
                method="POST",
                data=data,
                timeout_override=120,  # 120 second timeout for transcript queries
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

                        # Find MANE transcripts (both Select and Plus Clinical)
                        mane_select = None
                        mane_clinical = None
                        canonical_transcript_full = None
                        mane_select_full = None
                        mane_clinical_full = None

                        for transcript in transcripts:
                            transcript_id = transcript.get("id")

                            # Check if this is the canonical transcript
                            if transcript.get("is_canonical") == 1:
                                canonical_transcript_full = transcript

                            # Check for MANE transcripts
                            mane_list = transcript.get("MANE", [])
                            for mane_entry in mane_list:
                                mane_type = mane_entry.get("type")
                                if mane_type == "MANE_Select":
                                    mane_select = {
                                        "transcript_id": transcript_id,
                                        "refseq_match": mane_entry.get("refseq_match"),
                                        "transcript_start": transcript.get("start"),
                                        "transcript_end": transcript.get("end"),
                                        "biotype": transcript.get("biotype"),
                                        "length": transcript.get("length"),
                                    }
                                    mane_select_full = transcript
                                elif mane_type == "MANE_Plus_Clinical":
                                    mane_clinical = {
                                        "transcript_id": transcript_id,
                                        "refseq_match": mane_entry.get("refseq_match"),
                                        "transcript_start": transcript.get("start"),
                                        "transcript_end": transcript.get("end"),
                                        "biotype": transcript.get("biotype"),
                                        "length": transcript.get("length"),
                                    }
                                    mane_clinical_full = transcript

                        # Update the result with transcript data
                        gene_data = result[symbol]
                        if gene_data is not None:
                            gene_data["canonical_transcript"] = canonical_transcript
                            gene_data["mane_select"] = mane_select
                            gene_data["mane_clinical"] = mane_clinical
                            # Store full transcript data for coverage calculation
                            gene_data[
                                "canonical_transcript_full"
                            ] = canonical_transcript_full
                            gene_data["mane_select_full"] = mane_select_full
                            gene_data["mane_clinical_full"] = mane_clinical_full
                            # Store ALL transcripts with exon data for BED file generation
                            gene_data["all_transcripts"] = transcripts

        except requests.RequestException as e:
            logger.warning(f"Batch transcript request failed: {e}, using fallback")
            # Fallback to individual transcript requests for this batch
            for symbol, _gene_id in batch_items:
                gene_data = result[symbol]
                if gene_data is not None:
                    gene_data["canonical_transcript"] = None
                    gene_data["mane_select"] = None
                    gene_data["mane_clinical"] = None

    @functools.lru_cache(maxsize=1000)  # noqa: B019
    def rsid_to_coordinates(
        self,
        rsid: str,
        species: str = "human",
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
        self,
        gene_symbol: str,
        species: str = "human",
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
        self,
        gene_symbol: str,
        species: str = "human",
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
                "mane_select": None,
                "mane_clinical": None,
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
            "mane_select": gene_data.get("mane_select"),
            "mane_clinical": gene_data.get("mane_clinical"),
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

    def calculate_transcript_coverage(
        self,
        transcript_data: dict[str, Any],
        padding: int = 0,
    ) -> int | None:
        """
        Calculate transcript coverage including exons and padding.

        Args:
            transcript_data: Transcript data with exon information
            padding: Padding to add on both sides (in base pairs)

        Returns:
            Total coverage in base pairs or None if calculation fails
        """
        if not transcript_data or "Exon" not in transcript_data:
            return None

        try:
            exons = transcript_data["Exon"]
            total_exon_length = 0

            for exon in exons:
                start = exon.get("start")
                end = exon.get("end")
                if start is not None and end is not None:
                    total_exon_length += abs(end - start) + 1

            # Add padding on both sides
            total_coverage = total_exon_length + (2 * padding)
            return total_coverage

        except (KeyError, TypeError, ValueError):
            logger.warning("Failed to calculate transcript coverage for transcript")
            return None

    def calculate_gene_coverage(
        self,
        gene_start: int,
        gene_end: int,
        padding: int = 0,
    ) -> int | None:
        """
        Calculate gene coverage including padding.

        Args:
            gene_start: Gene start position
            gene_end: Gene end position
            padding: Padding to add on both sides (in base pairs)

        Returns:
            Total gene coverage in base pairs or None if calculation fails
        """
        if gene_start is None or gene_end is None:
            return None

        try:
            gene_length = abs(gene_end - gene_start) + 1
            total_coverage = gene_length + (2 * padding)
            return total_coverage

        except (TypeError, ValueError):
            logger.warning(
                f"Failed to calculate gene coverage for positions {gene_start}-{gene_end}",
            )
            return None

    def get_transcript_exons(
        self,
        transcript_id: str,
        species: str = "human",
    ) -> list[dict[str, Any]]:
        """
        Get exon coordinates for a specific transcript.

        Args:
            transcript_id: Ensembl transcript ID
            species: Species name (default: "human")

        Returns:
            List of exon dictionaries with coordinates
        """
        try:
            response = self._make_request(f"lookup/id/{transcript_id}?expand=1")

            if isinstance(response, dict) and "Exon" in response:
                exons = []
                for exon in response["Exon"]:
                    exon_info = {
                        "exon_id": exon.get("id"),
                        "chromosome": exon.get("seq_region_name"),
                        "start": exon.get("start"),
                        "end": exon.get("end"),
                        "strand": exon.get("strand"),
                        "rank": exon.get("rank", 0),
                    }
                    exons.append(exon_info)

                # Sort by rank to maintain exon order
                exons.sort(key=lambda x: x["rank"])
                return exons

        except (requests.RequestException, ValueError) as e:
            logger.warning(f"Failed to get exons for transcript {transcript_id}: {e}")

        return []

    def get_gene_exons_by_transcript_type(
        self,
        gene_data: dict[str, Any],
        transcript_type: str = "canonical",
    ) -> list[dict[str, Any]]:
        """
        Get exon coordinates for a gene using a specific transcript type.

        Args:
            gene_data: Gene data from batch annotation (with transcript info)
            transcript_type: Type of transcript ("canonical", "mane_select", "mane_clinical")

        Returns:
            List of exon dictionaries with coordinates and gene info
        """
        transcript_id = None

        # Get the appropriate transcript ID based on type
        if transcript_type == "canonical":
            transcript_info = gene_data.get("canonical_transcript")
            if transcript_info:
                transcript_id = transcript_info.get("transcript_id")
        elif transcript_type == "mane_select":
            mane_info = gene_data.get("mane_select")
            if mane_info:
                transcript_id = mane_info.get("transcript_id")
        elif transcript_type == "mane_clinical":
            mane_info = gene_data.get("mane_clinical")
            if mane_info:
                transcript_id = mane_info.get("transcript_id")

        if not transcript_id:
            return []

        # Get exons for this transcript
        exons = self.get_transcript_exons(transcript_id)

        # Add gene information to each exon
        gene_symbol = gene_data.get("gene_symbol", "")
        gene_id = gene_data.get("gene_id", "")

        for exon in exons:
            exon["gene_symbol"] = gene_symbol
            exon["gene_id"] = gene_id
            exon["transcript_id"] = transcript_id
            exon["transcript_type"] = transcript_type

        return exons

    def get_variations_batch(
        self,
        rsids: list[str],
        species: str = "homo_sapiens",
        batch_size: int = 25,
    ) -> dict[str, dict[str, Any]]:
        """
        Get variation information for multiple rsIDs using batched requests.

        This method automatically splits large requests into smaller batches
        to stay under the Ensembl API limit of 200 items per request.

        Args:
            rsids: List of rsID strings (e.g., ["rs123456", "rs789012"])
            species: Species name (default: "homo_sapiens")
            batch_size: Maximum number of rsIDs per batch (default: 25)

        Returns:
            Dictionary mapping rsID to variation information

        Raises:
            requests.RequestException: If request fails
        """
        if not rsids:
            return {}

        # Remove duplicates while preserving order
        unique_rsids = list(dict.fromkeys(rsids))

        # If we have fewer rsIDs than batch size, process in single batch
        if len(unique_rsids) <= batch_size:
            return self._get_variations_single_batch(unique_rsids, species)

        # Split into batches for large requests
        logger.info(
            f"🔍 Fetching variation data for {len(unique_rsids)} rsIDs from Ensembl in batches of {batch_size}",
            extra={
                "rsid_count": len(unique_rsids),
                "species": species,
                "batch_size": batch_size,
            },
        )

        all_results = {}
        batches = [
            unique_rsids[i : i + batch_size]
            for i in range(0, len(unique_rsids), batch_size)
        ]

        for batch_idx, batch_rsids in enumerate(batches):
            logger.debug(
                f"Processing batch {batch_idx + 1}/{len(batches)} with {len(batch_rsids)} rsIDs",
            )

            try:
                batch_results = self._get_variations_single_batch(batch_rsids, species)
                all_results.update(batch_results)

                logger.debug(
                    f"✅ Batch {batch_idx + 1}/{len(batches)} completed: {len(batch_results)} variations fetched",
                )

                # Add a small delay between batches to avoid overwhelming the API
                if batch_idx < len(batches) - 1:  # Don't delay after the last batch
                    time.sleep(0.5)  # 500ms delay between batches

            except requests.RequestException as e:
                error_msg = str(e)
                is_timeout = (
                    "timeout" in error_msg.lower() or "timed out" in error_msg.lower()
                )

                if is_timeout:
                    logger.warning(
                        f"⏱️ Batch {batch_idx + 1}/{len(batches)} timed out with {len(batch_rsids)} rsIDs. "
                        f"Consider reducing batch size further.",
                        extra={
                            "batch_rsids": batch_rsids,
                            "species": species,
                            "timeout_error": True,
                        },
                    )
                else:
                    logger.error(
                        f"❌ Batch {batch_idx + 1}/{len(batches)} failed: {error_msg}",
                        extra={
                            "batch_rsids": batch_rsids,
                            "species": species,
                        },
                    )
                # Continue with other batches even if one fails
                continue

        logger.info(
            f"✅ Successfully fetched {len(all_results)} variations from {len(batches)} batches",
            extra={
                "requested_count": len(unique_rsids),
                "returned_count": len(all_results),
                "batch_count": len(batches),
            },
        )

        return all_results

    def _get_variations_single_batch(
        self,
        rsids: list[str],
        species: str = "homo_sapiens",
    ) -> dict[str, dict[str, Any]]:
        """
        Get variation information for a single batch of rsIDs.

        Args:
            rsids: List of rsID strings (max 200)
            species: Species name (default: "homo_sapiens")

        Returns:
            Dictionary mapping rsID to variation information

        Raises:
            requests.RequestException: If request fails
        """
        if not rsids:
            return {}

        if len(rsids) > 200:
            logger.warning(
                f"Batch size {len(rsids)} exceeds Ensembl API limit of 200. "
                f"Consider using get_variations_batch() instead.",
            )

        url = f"variation/{species}/"
        data = {"ids": rsids}

        try:
            # Use longer timeout for batch variation requests as they take more time
            timeout_override = 90 if len(rsids) > 10 else None
            response = self._make_request(
                url,
                method="POST",
                data=data,
                timeout_override=timeout_override,
            )

            if isinstance(response, dict):
                return response
            else:
                logger.warning(f"Unexpected response format: {type(response)}")
                return {}

        except requests.RequestException as e:
            logger.error(
                f"❌ Failed to fetch variations from Ensembl: {e!s}",
                extra={
                    "rsid_count": len(rsids),
                    "species": species,
                },
            )
            raise

    def extract_coordinates_from_variation(
        self,
        variation_data: dict[str, Any],
        preferred_assembly: str = "GRCh38",
    ) -> dict[str, Any] | None:
        """
        Extract coordinate information from Ensembl variation data.

        Args:
            variation_data: Variation data from Ensembl API
            preferred_assembly: Preferred genome assembly (default: "GRCh38")

        Returns:
            Dictionary with coordinate information or None if not found
        """
        if not variation_data or "mappings" not in variation_data:
            return None

        mappings = variation_data["mappings"]
        if not mappings:
            return None

        # Find mapping for preferred assembly
        preferred_mapping = None
        fallback_mapping = None

        for mapping in mappings:
            if mapping.get("assembly_name") == preferred_assembly:
                preferred_mapping = mapping
                break
            elif mapping.get("coord_system") == "chromosome":
                fallback_mapping = mapping

        # Use preferred mapping or fallback to any chromosome mapping
        selected_mapping = preferred_mapping or fallback_mapping

        if not selected_mapping:
            return None

        # Extract coordinate information
        try:
            chromosome = selected_mapping["seq_region_name"]
            start = selected_mapping["start"]
            end = selected_mapping["end"]
            strand = selected_mapping.get("strand", 1)
            allele_string = selected_mapping.get("allele_string", "")
            assembly = selected_mapping.get("assembly_name", preferred_assembly)

            # Parse allele string to get ref and alt alleles
            ref_allele = ""
            alt_allele = ""
            if allele_string and "/" in allele_string:
                alleles = allele_string.split("/")
                if len(alleles) >= 2:
                    ref_allele = alleles[0]
                    alt_allele = alleles[1]  # Take first alternative allele

            return {
                "chromosome": str(chromosome),
                "start": int(start),
                "end": int(end),
                "strand": int(strand),
                "ref_allele": ref_allele,
                "alt_allele": alt_allele,
                "allele_string": allele_string,
                "assembly": assembly,
            }

        except (KeyError, ValueError, TypeError) as e:
            logger.debug(f"Failed to extract coordinates from mapping: {e}")
            return None

    def clear_cache(self) -> None:
        """Clear all cached results."""
        self.get_gene_coordinates.cache_clear()
        self.rsid_to_coordinates.cache_clear()
