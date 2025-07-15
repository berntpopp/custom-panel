"""
gnomAD API Client for variant resolution and coordinate liftover.

This module provides a comprehensive interface to the gnomAD GraphQL API
for resolving SNP identifiers and performing coordinate liftover between
genome builds.
"""

import json
import logging
import time
from typing import Any, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .gnomad_cache import GnomADCache

logger = logging.getLogger(__name__)


class PerformanceLogger:
    """Helper class for logging performance metrics."""

    def __init__(self, logger: logging.Logger):
        self.logger = logger

    def log_request_timing(
        self,
        operation: str,
        duration: float,
        success: bool,
        details: Optional[dict[str, Any]] = None,
    ) -> None:
        """Log API request timing with structured data."""
        status = "SUCCESS" if success else "FAILED"
        self.logger.info(
            f"🕐 {operation} completed in {duration:.3f}s [{status}]",
            extra={
                "operation": operation,
                "duration_seconds": duration,
                "success": success,
                "details": details or {},
            },
        )

    def log_cache_performance(
        self, operation: str, cache_hit: bool, details: Optional[dict[str, Any]] = None
    ) -> None:
        """Log cache performance metrics."""
        hit_status = "HIT" if cache_hit else "MISS"
        emoji = "💾" if cache_hit else "🔍"
        self.logger.debug(
            f"{emoji} Cache {hit_status} for {operation}",
            extra={
                "operation": operation,
                "cache_hit": cache_hit,
                "details": details or {},
            },
        )

    def log_error_with_context(
        self, operation: str, error: Exception, context: Optional[dict[str, Any]] = None
    ) -> None:
        """Log errors with structured context."""
        self.logger.error(
            f"❌ {operation} failed: {error}",
            extra={
                "operation": operation,
                "error_type": type(error).__name__,
                "error_message": str(error),
                "context": context or {},
            },
            exc_info=True,
        )


class RateLimiter:
    """Rate limiter for API calls with exponential backoff support."""

    def __init__(self, requests_per_second: float = 5.0):
        self.requests_per_second = requests_per_second
        self.last_request_time = 0.0
        self.min_interval = 1.0 / requests_per_second
        self.consecutive_429_errors = 0
        self.backoff_multiplier = 2.0
        self.max_backoff = 60.0  # Maximum backoff time in seconds
        self.base_429_delay = 1.0  # Base delay for 429 errors

    def wait_if_needed(self) -> None:
        """Wait if needed to respect rate limit."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time

        # Calculate minimum interval including exponential backoff for 429 errors
        effective_interval = self.min_interval
        if self.consecutive_429_errors > 0:
            backoff_delay = min(
                self.base_429_delay
                * (self.backoff_multiplier**self.consecutive_429_errors),
                self.max_backoff,
            )
            effective_interval = max(self.min_interval, backoff_delay)
            logger.info(
                f"🐌 Exponential backoff active due to {self.consecutive_429_errors} consecutive 429 errors, waiting {effective_interval:.2f}s",
                extra={
                    "consecutive_429_errors": self.consecutive_429_errors,
                    "backoff_delay": backoff_delay,
                    "effective_interval": effective_interval,
                },
            )

        if time_since_last < effective_interval:
            sleep_time = effective_interval - time_since_last
            logger.debug(f"Rate limiting: sleeping {sleep_time:.2f}s")
            time.sleep(sleep_time)

        self.last_request_time = time.time()

    def record_429_error(self) -> None:
        """Record a 429 error for exponential backoff calculation."""
        self.consecutive_429_errors += 1
        logger.warning(
            f"⚠️ 429 Rate limit error #{self.consecutive_429_errors}, increasing backoff",
            extra={"consecutive_429_errors": self.consecutive_429_errors},
        )

    def reset_429_errors(self) -> None:
        """Reset 429 error counter after successful request."""
        if self.consecutive_429_errors > 0:
            logger.info(
                f"✅ Request successful, resetting 429 error counter from {self.consecutive_429_errors}",
                extra={"previous_error_count": self.consecutive_429_errors},
            )
            self.consecutive_429_errors = 0


class GnomADAPIError(Exception):
    """Custom exception for gnomAD API errors."""

    pass


class GnomADClient:
    """
    Client for interacting with the gnomAD GraphQL API.

    Provides methods for:
    - Resolving rsIDs from coordinates
    - Performing coordinate liftover between genome builds
    - Validating variants against gnomAD database
    """

    def __init__(
        self,
        base_url: str = "https://gnomad.broadinstitute.org/api",
        rate_limit: float = 5.0,
        timeout: int = 30,
        retry_attempts: int = 3,
        cache_dir: str = ".cache/gnomad",
        cache_ttl_days: int = 30,
    ):
        """
        Initialize gnomAD client.

        Args:
            base_url: Base URL for gnomAD API
            rate_limit: Requests per second limit
            timeout: Request timeout in seconds
            retry_attempts: Number of retry attempts
            cache_dir: Directory for caching responses
            cache_ttl_days: Cache TTL in days
        """
        self.base_url = base_url
        self.timeout = timeout
        self.retry_attempts = retry_attempts

        # Initialize rate limiter
        self.rate_limiter = RateLimiter(rate_limit)

        # Initialize cache
        self.cache = GnomADCache(cache_dir, cache_ttl_days)

        # Setup HTTP session with retries
        self.session = requests.Session()

        # Custom retry strategy with exponential backoff for 429 errors
        retry_strategy = Retry(
            total=retry_attempts,
            backoff_factor=1.5,  # Increased backoff factor for better handling
            status_forcelist=[500, 502, 503, 504],  # Removed 429 to handle manually
            allowed_methods=["GET", "POST"],
            raise_on_status=False,  # Don't raise immediately, let us handle 429s
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

        # Set headers
        self.session.headers.update(
            {
                "Content-Type": "application/json",
                "User-Agent": "custom-panel-gnomad-client/1.0",
            }
        )

        # Initialize performance logger
        self.perf_logger = PerformanceLogger(logger)

        logger.info(
            f"🚀 Initialized gnomAD client with base URL: {base_url}",
            extra={
                "base_url": base_url,
                "rate_limit": rate_limit,
                "timeout": timeout,
                "retry_attempts": retry_attempts,
                "cache_dir": cache_dir,
                "cache_ttl_days": cache_ttl_days,
            },
        )

    def _make_graphql_request(
        self, query: str, variables: dict[str, Any], max_429_retries: int = 5
    ) -> dict[str, Any]:
        """
        Make a GraphQL request to gnomAD API with exponential backoff for 429 errors.

        Args:
            query: GraphQL query string
            variables: Query variables
            max_429_retries: Maximum number of retries for 429 errors

        Returns:
            Response data

        Raises:
            GnomADAPIError: If request fails
        """
        operation = f"GraphQL {query.split()[1].split('(')[0]}"

        for attempt in range(max_429_retries + 1):
            start_time = time.time()

            # Apply rate limiting (includes exponential backoff for 429s)
            self.rate_limiter.wait_if_needed()

            payload = {
                "query": query,
                "variables": variables,
            }

            try:
                response = self.session.post(
                    self.base_url,
                    json=payload,
                    timeout=self.timeout,
                )

                # Handle 429 errors specifically with exponential backoff
                if response.status_code == 429:
                    duration = time.time() - start_time
                    if attempt < max_429_retries:
                        self.rate_limiter.record_429_error()
                        self.perf_logger.log_request_timing(
                            operation,
                            duration,
                            False,
                            {
                                "status_code": 429,
                                "attempt": attempt + 1,
                                "max_retries": max_429_retries,
                                "variables": variables,
                                "will_retry": True,
                            },
                        )
                        logger.warning(
                            f"🔄 429 Rate limit hit, retrying {attempt + 1}/{max_429_retries} for {operation}",
                            extra={
                                "attempt": attempt + 1,
                                "max_retries": max_429_retries,
                                "operation": operation,
                                "variables": variables,
                            },
                        )
                        continue  # Retry with exponential backoff
                    else:
                        # Final attempt failed with 429
                        self.perf_logger.log_request_timing(
                            operation,
                            duration,
                            False,
                            {
                                "status_code": 429,
                                "attempt": attempt + 1,
                                "max_retries": max_429_retries,
                                "variables": variables,
                                "will_retry": False,
                            },
                        )
                        raise GnomADAPIError(
                            f"Request failed: {response.status_code} {response.reason} "
                            f"(too many 429 error responses after {max_429_retries} retries)"
                        )

                # For non-429 status codes, use normal handling
                response.raise_for_status()

                data = response.json()
                duration = time.time() - start_time

                # Reset 429 error counter on successful request
                self.rate_limiter.reset_429_errors()

                # Check for GraphQL errors
                if "errors" in data:
                    error_messages = [
                        error.get("message", "Unknown error")
                        for error in data["errors"]
                    ]
                    self.perf_logger.log_request_timing(
                        operation,
                        duration,
                        False,
                        {"error_count": len(data["errors"]), "variables": variables},
                    )
                    raise GnomADAPIError(f"GraphQL errors: {'; '.join(error_messages)}")

                self.perf_logger.log_request_timing(
                    operation,
                    duration,
                    True,
                    {"response_size": len(str(data)), "variables": variables},
                )
                return data.get("data", {})

            except requests.exceptions.RequestException as e:
                duration = time.time() - start_time

                # Check if it's a connection pool error indicating too many 429s
                if "too many 429 error responses" in str(e):
                    self.perf_logger.log_request_timing(
                        operation, duration, False, {"variables": variables}
                    )
                    self.perf_logger.log_error_with_context(
                        operation,
                        e,
                        {
                            "variables": variables,
                            "url": self.base_url,
                            "attempt": attempt + 1,
                        },
                    )
                    if attempt < max_429_retries:
                        self.rate_limiter.record_429_error()
                        logger.warning(
                            f"🔄 Connection pool 429 error, retrying {attempt + 1}/{max_429_retries} for {operation}",
                            extra={
                                "attempt": attempt + 1,
                                "max_retries": max_429_retries,
                                "operation": operation,
                                "error": str(e),
                            },
                        )
                        continue  # Retry with exponential backoff
                    else:
                        raise GnomADAPIError(f"Request failed: {str(e)}") from e
                else:
                    # Other request exceptions, don't retry
                    self.perf_logger.log_request_timing(
                        operation, duration, False, {"variables": variables}
                    )
                    self.perf_logger.log_error_with_context(
                        operation, e, {"variables": variables, "url": self.base_url}
                    )
                    raise GnomADAPIError(f"Request failed: {str(e)}") from e

            except json.JSONDecodeError as e:
                duration = time.time() - start_time
                self.perf_logger.log_request_timing(
                    operation, duration, False, {"variables": variables}
                )
                self.perf_logger.log_error_with_context(
                    operation,
                    e,
                    {
                        "variables": variables,
                        "response_text": getattr(response, "text", "N/A")[:500],
                    },
                )
                raise GnomADAPIError(f"Invalid JSON response: {str(e)}") from e

        # This should never be reached, but just in case
        raise GnomADAPIError(
            f"Max retries ({max_429_retries}) exceeded for {operation}"
        )

    def resolve_rsid(
        self,
        chromosome: str,
        position: int,
        ref_allele: str,
        alt_allele: str,
        build: str = "GRCh38",
    ) -> Optional[str]:
        """
        Resolve rsID from genomic coordinates.

        Args:
            chromosome: Chromosome (e.g., "1", "X")
            position: Genomic position
            ref_allele: Reference allele
            alt_allele: Alternate allele
            build: Genome build ("GRCh37" or "GRCh38")

        Returns:
            rsID if found, None otherwise
        """
        # Create variant ID
        variant_id = f"{chromosome}-{position}-{ref_allele}-{alt_allele}"

        # Check cache first
        cache_key = self.cache.get_rsid_cache_key(
            chromosome, position, ref_allele, alt_allele, build
        )
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            self.perf_logger.log_cache_performance(
                f"rsID resolution: {variant_id}", True, {"build": build}
            )
            return cached_result

        self.perf_logger.log_cache_performance(
            f"rsID resolution: {variant_id}", False, {"build": build}
        )

        # Determine dataset based on build
        dataset = "gnomad_r2_1" if build == "GRCh37" else "gnomad_r4"

        query = """
        query variant($variantId: String!, $dataset: DatasetId!) {
            variant(variantId: $variantId, dataset: $dataset) {
                variant_id
                rsid
                chrom
                pos
                ref
                alt
            }
        }
        """

        variables = {
            "variantId": variant_id,
            "dataset": dataset,
        }

        try:
            logger.debug(
                f"🔍 Resolving rsID for variant: {variant_id} (build: {build})",
                extra={"variant_id": variant_id, "build": build, "dataset": dataset},
            )
            data = self._make_graphql_request(query, variables)

            variant = data.get("variant")
            if variant:
                rsid = variant.get("rsid")
                if rsid:
                    logger.debug(
                        f"✅ Resolved rsID: {variant_id} -> {rsid}",
                        extra={"variant_id": variant_id, "rsid": rsid, "build": build},
                    )
                    self.cache.set(cache_key, rsid)
                    return rsid
                else:
                    logger.debug(
                        f"❌ No rsID found for variant: {variant_id}",
                        extra={"variant_id": variant_id, "build": build},
                    )
                    self.cache.set(cache_key, None)
                    return None
            else:
                logger.debug(
                    f"❌ Variant not found in gnomAD: {variant_id}",
                    extra={
                        "variant_id": variant_id,
                        "build": build,
                        "dataset": dataset,
                    },
                )
                self.cache.set(cache_key, None)
                return None

        except GnomADAPIError as e:
            self.perf_logger.log_error_with_context(
                f"rsID resolution: {variant_id}",
                e,
                {"variant_id": variant_id, "build": build, "dataset": dataset},
            )
            return None

    def liftover_coordinates(
        self,
        source_variant_id: str,
        source_build: str,
        target_build: str,
    ) -> Optional[tuple[str, str, int, str, str]]:
        """
        Perform coordinate liftover between genome builds.

        Args:
            source_variant_id: Source variant ID (chr-pos-ref-alt)
            source_build: Source genome build ("GRCh37" or "GRCh38")
            target_build: Target genome build ("GRCh37" or "GRCh38")

        Returns:
            Tuple of (chromosome, target_variant_id, position, ref, alt) if successful, None otherwise
        """
        if source_build == target_build:
            logger.debug(
                f"Same build requested, returning original: {source_variant_id}"
            )
            # Parse the variant ID to return consistent format
            parts = source_variant_id.split("-")
            if len(parts) == 4:
                return parts[0], source_variant_id, int(parts[1]), parts[2], parts[3]
            return None

        # Check cache first
        cache_key = self.cache.get_liftover_cache_key(
            source_variant_id, source_build, target_build
        )
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            logger.debug(
                f"Cache hit for liftover: {source_variant_id} ({source_build} -> {target_build})"
            )
            return cached_result

        # Determine reference genome parameter
        reference_genome = "GRCh37" if source_build == "GRCh37" else "GRCh38"

        query = """
        query liftover($source_variant_id: String!, $reference_genome: ReferenceGenomeId!) {
            liftover(source_variant_id: $source_variant_id, reference_genome: $reference_genome) {
                source {
                    variant_id
                    reference_genome
                }
                liftover {
                    variant_id
                    reference_genome
                }
                datasets
            }
        }
        """

        variables = {
            "source_variant_id": source_variant_id,
            "reference_genome": reference_genome,
        }

        try:
            logger.debug(
                f"Lifting over coordinates: {source_variant_id} ({source_build} -> {target_build})"
            )
            data = self._make_graphql_request(query, variables)

            liftover_data = data.get("liftover")
            if liftover_data and liftover_data.get("liftover"):
                target_variant_id = liftover_data["liftover"]["variant_id"]

                # Parse target variant ID
                parts = target_variant_id.split("-")
                if len(parts) == 4:
                    chromosome, position, ref, alt = parts
                    result = (chromosome, target_variant_id, int(position), ref, alt)
                    logger.debug(
                        f"Liftover successful: {source_variant_id} -> {target_variant_id}"
                    )
                    self.cache.set(cache_key, result)
                    return result
                else:
                    logger.warning(
                        f"Invalid target variant ID format: {target_variant_id}"
                    )
                    self.cache.set(cache_key, None)
                    return None
            else:
                logger.debug(f"Liftover failed for variant: {source_variant_id}")
                self.cache.set(cache_key, None)
                return None

        except GnomADAPIError as e:
            logger.error(f"Failed to liftover coordinates for {source_variant_id}: {e}")
            return None

    def validate_variant(
        self,
        chromosome: str,
        position: int,
        ref_allele: str,
        alt_allele: str,
        build: str = "GRCh38",
    ) -> bool:
        """
        Validate that a variant exists in gnomAD.

        Args:
            chromosome: Chromosome
            position: Genomic position
            ref_allele: Reference allele
            alt_allele: Alternate allele
            build: Genome build

        Returns:
            True if variant exists, False otherwise
        """
        variant_id = f"{chromosome}-{position}-{ref_allele}-{alt_allele}"

        # Check cache first
        cache_key = f"validate:{build}:{variant_id}"
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            return cached_result

        # Determine dataset based on build
        dataset = "gnomad_r2_1" if build == "GRCh37" else "gnomad_r4"

        query = """
        query variant($variantId: String!, $dataset: DatasetId!) {
            variant(variantId: $variantId, dataset: $dataset) {
                variant_id
                chrom
                pos
                ref
                alt
            }
        }
        """

        variables = {
            "variantId": variant_id,
            "dataset": dataset,
        }

        try:
            data = self._make_graphql_request(query, variables)
            exists = data.get("variant") is not None
            self.cache.set(cache_key, exists)
            return exists

        except GnomADAPIError as e:
            logger.error(f"Failed to validate variant {variant_id}: {e}")
            return False

    def get_variant_by_rsid(self, rsid: str, build: str) -> Optional[dict[str, Any]]:
        """
        Get variant information by rsID.

        Args:
            rsid: rsID to look up
            build: Genome build

        Returns:
            Variant information dictionary or None if not found
        """
        # Check cache first
        cache_key = f"rsid_variant:{build}:{rsid}"
        cached_result = self.cache.get(cache_key)
        if cached_result is not None:
            return cached_result

        # Determine dataset based on build
        dataset = "gnomad_r2_1" if build == "GRCh37" else "gnomad_r4"

        query = """
        query variantByRsid($rsid: String!, $dataset: DatasetId!) {
            variants(rsid: $rsid, dataset: $dataset) {
                variant_id
                rsid
                chrom
                pos
                ref
                alt
            }
        }
        """

        variables = {
            "rsid": rsid,
            "dataset": dataset,
        }

        try:
            data = self._make_graphql_request(query, variables)
            variants = data.get("variants", [])

            # If we have multiple variants, take the first one (they're at the same position)
            if variants:
                variant_info = variants[0]
                if len(variants) > 1:
                    logger.debug(
                        f"Multiple variants found for rsID {rsid}, using first variant: {variant_info.get('variant_id')}",
                        extra={
                            "rsid": rsid,
                            "variant_count": len(variants),
                            "selected_variant": variant_info.get("variant_id"),
                        },
                    )
                self.cache.set(cache_key, variant_info)
                return variant_info
            else:
                # No variants found
                self.cache.set(cache_key, None)
                return None

        except GnomADAPIError as e:
            logger.error(f"Failed to get variant by rsID {rsid}: {e}")
            return None

    def get_cache_stats(self) -> dict[str, Any]:
        """Get cache statistics."""
        return self.cache.get_stats()

    def clear_cache(self) -> int:
        """Clear all cached data."""
        return self.cache.clear()
