"""
Cache manager for API responses.

This module provides a file-based caching system for API responses to improve
performance on repeated runs.
"""

import hashlib
import json
import logging
import time
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class CacheManager:
    """Manages file-based caching of API responses."""

    def __init__(
        self,
        cache_dir: str = ".cache",
        cache_ttl: int = 31536000,  # 365 days in seconds
        enabled: bool = True,
    ):
        """
        Initialize the cache manager.

        Args:
            cache_dir: Directory to store cache files
            cache_ttl: Time-to-live for cache entries in seconds
            enabled: Whether caching is enabled
        """
        self.cache_dir = Path(cache_dir)
        self.cache_ttl = cache_ttl
        self.enabled = enabled

        if self.enabled:
            self._ensure_cache_dir()

    def _ensure_cache_dir(self) -> None:
        """Ensure the cache directory exists."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories for different cache types
        (self.cache_dir / "ensembl").mkdir(exist_ok=True)
        (self.cache_dir / "hgnc").mkdir(exist_ok=True)
        (self.cache_dir / "clingen").mkdir(exist_ok=True)
        (self.cache_dir / "gencc").mkdir(exist_ok=True)
        (self.cache_dir / "acmg").mkdir(exist_ok=True)
        (self.cache_dir / "panelapp").mkdir(exist_ok=True)
        (self.cache_dir / "omim").mkdir(exist_ok=True)
        (self.cache_dir / "hpo").mkdir(exist_ok=True)
        (self.cache_dir / "cosmic").mkdir(exist_ok=True)

    def _generate_cache_key(
        self,
        service: str,
        endpoint: str,
        method: str = "GET",
        data: dict[str, Any] | None = None,
    ) -> str:
        """
        Generate a unique cache key for a request.

        Args:
            service: Service name (e.g., "ensembl", "hgnc")
            endpoint: API endpoint
            method: HTTP method
            data: Request data

        Returns:
            Cache key string
        """
        # Create a unique string from request parameters
        key_parts = [service, endpoint, method]

        if data:
            # Sort data keys for consistent hashing
            sorted_data = json.dumps(data, sort_keys=True)
            key_parts.append(sorted_data)

        key_string = "|".join(key_parts)

        # Generate SHA256 hash for the key
        key_hash = hashlib.sha256(key_string.encode()).hexdigest()

        # Use first 16 chars of hash + readable suffix
        readable_part = endpoint.replace("/", "_").replace("?", "_")[-20:]
        return f"{key_hash[:16]}_{readable_part}"

    def _get_cache_path(self, service: str, cache_key: str) -> Path:
        """Get the full path for a cache file."""
        return self.cache_dir / service / f"{cache_key}.json"

    def get(
        self,
        service: str,
        endpoint: str,
        method: str = "GET",
        data: dict[str, Any] | None = None,
    ) -> Any | None:
        """
        Retrieve cached response if available and valid.

        Args:
            service: Service name
            endpoint: API endpoint
            method: HTTP method
            data: Request data

        Returns:
            Cached response or None if not found/expired
        """
        if not self.enabled:
            return None

        cache_key = self._generate_cache_key(service, endpoint, method, data)
        cache_path = self._get_cache_path(service, cache_key)

        if not cache_path.exists():
            logger.debug(f"Cache miss for {service} {endpoint}")
            return None

        try:
            with open(cache_path) as f:
                cache_entry = json.load(f)

            # Check if cache is expired
            if time.time() - cache_entry["timestamp"] > self.cache_ttl:
                logger.info(f"⏰ Cache expired for {service} {endpoint}")
                cache_path.unlink()  # Remove expired cache
                return None

            # Show cache hits at INFO level so users can see the benefit
            if data and "symbols" in data:
                gene_count = len(data["symbols"])
                logger.info(f"🚀 Cache hit: {service} {endpoint} ({gene_count} genes)")
            elif data and "ids" in data:
                gene_count = len(data["ids"])
                logger.info(
                    f"🚀 Cache hit: {service} {endpoint} ({gene_count} gene IDs)"
                )
            else:
                logger.info(f"🚀 Cache hit: {service} {endpoint}")
            return cache_entry["response"]

        except (OSError, json.JSONDecodeError, KeyError) as e:
            logger.warning(f"Error reading cache for {cache_key}: {e}")
            # Remove corrupted cache file
            if cache_path.exists():
                cache_path.unlink()
            return None

    def set(
        self,
        service: str,
        endpoint: str,
        method: str,
        data: dict[str, Any] | None,
        response: Any,
    ) -> None:
        """
        Cache a response.

        Args:
            service: Service name
            endpoint: API endpoint
            method: HTTP method
            data: Request data
            response: Response to cache
        """
        if not self.enabled:
            return

        cache_key = self._generate_cache_key(service, endpoint, method, data)
        cache_path = self._get_cache_path(service, cache_key)

        cache_entry = {
            "timestamp": time.time(),
            "endpoint": endpoint,
            "method": method,
            "response": response,
        }

        try:
            # Write atomically to avoid corruption
            temp_path = cache_path.with_suffix(".tmp")
            with open(temp_path, "w") as f:
                json.dump(cache_entry, f, indent=2)

            # Atomic rename
            temp_path.replace(cache_path)

            # Show cache saves at INFO level for transparency
            if data and "symbols" in data:
                gene_count = len(data["symbols"])
                logger.info(
                    f"💾 Cached response: {service} {endpoint} ({gene_count} genes)"
                )
            elif data and "ids" in data:
                gene_count = len(data["ids"])
                logger.info(
                    f"💾 Cached response: {service} {endpoint} ({gene_count} gene IDs)"
                )
            else:
                logger.info(f"💾 Cached response: {service} {endpoint}")

        except OSError as e:
            logger.warning(f"Failed to cache response for {cache_key}: {e}")
            # Clean up temp file if it exists
            if temp_path.exists():
                temp_path.unlink()

    def clear(self, service: str | None = None) -> int:
        """
        Clear cache entries.

        Args:
            service: Service to clear (None for all)

        Returns:
            Number of entries cleared
        """
        count = 0

        if service:
            service_dir = self.cache_dir / service
            if service_dir.exists():
                for cache_file in service_dir.glob("*.json"):
                    cache_file.unlink()
                    count += 1
        else:
            # Clear all services
            for service_dir in self.cache_dir.iterdir():
                if service_dir.is_dir():
                    for cache_file in service_dir.glob("*.json"):
                        cache_file.unlink()
                        count += 1

        logger.info(f"Cleared {count} cache entries")
        return count

    def get_stats(self) -> dict[str, Any]:
        """
        Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        stats: dict[str, Any] = {
            "enabled": self.enabled,
            "cache_dir": str(self.cache_dir),
            "cache_ttl_days": self.cache_ttl / 86400,
            "services": {},
        }

        if self.enabled and self.cache_dir.exists():
            total_size = 0
            total_files = 0

            for service_dir in self.cache_dir.iterdir():
                if service_dir.is_dir():
                    service_files = list(service_dir.glob("*.json"))
                    service_size = sum(f.stat().st_size for f in service_files)

                    stats["services"][service_dir.name] = {
                        "files": len(service_files),
                        "size_mb": round(service_size / 1024 / 1024, 2),
                    }

                    total_files += len(service_files)
                    total_size += service_size

            stats["total_files"] = total_files
            stats["total_size_mb"] = round(total_size / 1024 / 1024, 2)

        return stats

    def cleanup_expired(self) -> int:
        """
        Remove expired cache entries.

        Returns:
            Number of entries removed
        """
        if not self.enabled:
            return 0

        count = 0
        current_time = time.time()

        for service_dir in self.cache_dir.iterdir():
            if service_dir.is_dir():
                for cache_file in service_dir.glob("*.json"):
                    try:
                        with open(cache_file) as f:
                            cache_entry = json.load(f)

                        if current_time - cache_entry["timestamp"] > self.cache_ttl:
                            cache_file.unlink()
                            count += 1

                    except (OSError, json.JSONDecodeError, KeyError):
                        # Remove corrupted files
                        cache_file.unlink()
                        count += 1

        if count > 0:
            logger.info(f"Cleaned up {count} expired cache entries")

        return count
