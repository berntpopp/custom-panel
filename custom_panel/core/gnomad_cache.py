"""
Cache management for gnomAD API responses.

This module provides a comprehensive caching system for gnomAD API responses
with TTL-based expiration, persistent storage, and cache statistics.
"""

import json
import logging
import pickle
import time
from datetime import datetime
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class GnomADCache:
    """
    Persistent cache for gnomAD API responses with TTL-based expiration.

    Provides separate caches for different operation types:
    - rsID resolution
    - Coordinate liftover
    - Variant validation
    """

    def __init__(self, cache_dir: str = ".cache/gnomad", ttl_days: int = 30):
        """
        Initialize cache manager.

        Args:
            cache_dir: Directory for cache files
            ttl_days: Time-to-live in days
        """
        self.cache_dir = Path(cache_dir)
        self.ttl_seconds = ttl_days * 24 * 60 * 60  # Convert to seconds

        # Create cache directory if it doesn't exist
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Cache file paths
        self.rsid_cache_file = self.cache_dir / "rsid_cache.pkl"
        self.liftover_cache_file = self.cache_dir / "liftover_cache.pkl"
        self.validation_cache_file = self.cache_dir / "validation_cache.pkl"
        self.stats_file = self.cache_dir / "cache_stats.json"

        # Load existing caches
        self.rsid_cache = self._load_cache(self.rsid_cache_file)
        self.liftover_cache = self._load_cache(self.liftover_cache_file)
        self.validation_cache = self._load_cache(self.validation_cache_file)

        # Initialize stats
        self.stats = self._load_stats()

        logger.info(f"Initialized gnomAD cache in {cache_dir} with TTL {ttl_days} days")

    def _load_cache(self, cache_file: Path) -> dict[str, tuple[Any, float]]:
        """
        Load cache from file.

        Args:
            cache_file: Path to cache file

        Returns:
            Cache dictionary with (value, timestamp) tuples
        """
        if cache_file.exists():
            try:
                with open(cache_file, "rb") as f:
                    cache = pickle.load(f)
                logger.debug(
                    f"Loaded cache from {cache_file} with {len(cache)} entries"
                )
                return cache
            except Exception as e:
                logger.warning(f"Failed to load cache from {cache_file}: {e}")
                return {}
        return {}

    def _save_cache(
        self, cache: dict[str, tuple[Any, float]], cache_file: Path
    ) -> None:
        """
        Save cache to file.

        Args:
            cache: Cache dictionary
            cache_file: Path to cache file
        """
        try:
            with open(cache_file, "wb") as f:
                pickle.dump(cache, f)
            logger.debug(f"Saved cache to {cache_file} with {len(cache)} entries")
        except Exception as e:
            logger.error(f"Failed to save cache to {cache_file}: {e}")

    def _load_stats(self) -> dict[str, Any]:
        """Load cache statistics."""
        if self.stats_file.exists():
            try:
                with open(self.stats_file) as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load cache stats: {e}")

        return {
            "hits": 0,
            "misses": 0,
            "sets": 0,
            "expired": 0,
            "last_cleanup": None,
        }

    def _save_stats(self) -> None:
        """Save cache statistics."""
        try:
            with open(self.stats_file, "w") as f:
                json.dump(self.stats, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save cache stats: {e}")

    def _is_expired(self, timestamp: float) -> bool:
        """Check if cache entry is expired."""
        return time.time() - timestamp > self.ttl_seconds

    def _cleanup_expired(self, cache: dict[str, tuple[Any, float]]) -> int:
        """Remove expired entries from cache."""
        expired_keys = [
            key for key, (_, timestamp) in cache.items() if self._is_expired(timestamp)
        ]

        for key in expired_keys:
            del cache[key]

        if expired_keys:
            self.stats["expired"] += len(expired_keys)
            logger.debug(f"Cleaned up {len(expired_keys)} expired cache entries")

        return len(expired_keys)

    def get_rsid_cache_key(
        self,
        chromosome: str,
        position: int,
        ref_allele: str,
        alt_allele: str,
        build: str,
    ) -> str:
        """Generate cache key for rsID resolution."""
        return f"rsid:{build}:{chromosome}:{position}:{ref_allele}:{alt_allele}"

    def get_liftover_cache_key(
        self,
        variant_id: str,
        source_build: str,
        target_build: str,
    ) -> str:
        """Generate cache key for coordinate liftover."""
        return f"liftover:{source_build}:{target_build}:{variant_id}"

    def get(self, key: str) -> Any:
        """
        Get value from cache.

        Args:
            key: Cache key

        Returns:
            Cached value or None if not found/expired
        """
        # Determine which cache to use based on key prefix
        if key.startswith("rsid:") or key.startswith("rsid_variant:"):
            cache = self.rsid_cache
            cache_file = self.rsid_cache_file
        elif key.startswith("liftover:"):
            cache = self.liftover_cache
            cache_file = self.liftover_cache_file
        elif key.startswith("validate:"):
            cache = self.validation_cache
            cache_file = self.validation_cache_file
        else:
            logger.warning(f"Unknown cache key prefix: {key}")
            return None

        if key in cache:
            value, timestamp = cache[key]
            if self._is_expired(timestamp):
                del cache[key]
                self.stats["expired"] += 1
                self._save_cache(cache, cache_file)
                self.stats["misses"] += 1
                return None
            else:
                self.stats["hits"] += 1
                return value

        self.stats["misses"] += 1
        return None

    def set(self, key: str, value: Any) -> None:
        """
        Set value in cache.

        Args:
            key: Cache key
            value: Value to cache
        """
        timestamp = time.time()

        # Determine which cache to use based on key prefix
        if key.startswith("rsid:") or key.startswith("rsid_variant:"):
            cache = self.rsid_cache
            cache_file = self.rsid_cache_file
        elif key.startswith("liftover:"):
            cache = self.liftover_cache
            cache_file = self.liftover_cache_file
        elif key.startswith("validate:"):
            cache = self.validation_cache
            cache_file = self.validation_cache_file
        else:
            logger.warning(f"Unknown cache key prefix: {key}")
            return

        cache[key] = (value, timestamp)
        self.stats["sets"] += 1

        # Save cache to disk
        self._save_cache(cache, cache_file)
        self._save_stats()

    def cleanup_expired(self) -> int:
        """
        Clean up expired entries from all caches.

        Returns:
            Number of expired entries removed
        """
        total_expired = 0

        # Clean up each cache
        total_expired += self._cleanup_expired(self.rsid_cache)
        total_expired += self._cleanup_expired(self.liftover_cache)
        total_expired += self._cleanup_expired(self.validation_cache)

        # Save cleaned caches
        if total_expired > 0:
            self._save_cache(self.rsid_cache, self.rsid_cache_file)
            self._save_cache(self.liftover_cache, self.liftover_cache_file)
            self._save_cache(self.validation_cache, self.validation_cache_file)

        # Update stats
        self.stats["last_cleanup"] = datetime.now().isoformat()
        self._save_stats()

        return total_expired

    def clear(self) -> int:
        """
        Clear all cached data.

        Returns:
            Number of entries removed
        """
        total_entries = (
            len(self.rsid_cache) + len(self.liftover_cache) + len(self.validation_cache)
        )

        # Clear all caches
        self.rsid_cache.clear()
        self.liftover_cache.clear()
        self.validation_cache.clear()

        # Save empty caches
        self._save_cache(self.rsid_cache, self.rsid_cache_file)
        self._save_cache(self.liftover_cache, self.liftover_cache_file)
        self._save_cache(self.validation_cache, self.validation_cache_file)

        # Reset stats
        self.stats = {
            "hits": 0,
            "misses": 0,
            "sets": 0,
            "expired": 0,
            "last_cleanup": datetime.now().isoformat(),
        }
        self._save_stats()

        logger.info(f"Cleared cache, removed {total_entries} entries")
        return total_entries

    def get_stats(self) -> dict[str, Any]:
        """
        Get cache statistics.

        Returns:
            Dictionary with cache statistics
        """
        # Calculate hit rate
        total_requests = self.stats["hits"] + self.stats["misses"]
        hit_rate = (
            (self.stats["hits"] / total_requests * 100) if total_requests > 0 else 0
        )

        # Get cache sizes
        cache_sizes = {
            "rsid_cache": len(self.rsid_cache),
            "liftover_cache": len(self.liftover_cache),
            "validation_cache": len(self.validation_cache),
        }

        # Get cache file sizes
        cache_file_sizes = {}
        for name, file_path in [
            ("rsid_cache", self.rsid_cache_file),
            ("liftover_cache", self.liftover_cache_file),
            ("validation_cache", self.validation_cache_file),
        ]:
            if file_path.exists():
                cache_file_sizes[name] = file_path.stat().st_size
            else:
                cache_file_sizes[name] = 0

        return {
            **self.stats,
            "hit_rate": round(hit_rate, 2),
            "total_requests": total_requests,
            "cache_sizes": cache_sizes,
            "total_cached_entries": sum(cache_sizes.values()),
            "cache_file_sizes_bytes": cache_file_sizes,
            "total_cache_size_bytes": sum(cache_file_sizes.values()),
            "ttl_days": self.ttl_seconds / (24 * 60 * 60),
        }

    def get_cache_info(self) -> str:
        """Get formatted cache information."""
        stats = self.get_stats()

        info = f"""
gnomAD Cache Information:
-------------------------
Total Requests: {stats['total_requests']}
Cache Hits: {stats['hits']} ({stats['hit_rate']}%)
Cache Misses: {stats['misses']}
Cache Sets: {stats['sets']}
Expired Entries: {stats['expired']}

Cache Sizes:
- rsID Cache: {stats['cache_sizes']['rsid_cache']} entries
- Liftover Cache: {stats['cache_sizes']['liftover_cache']} entries
- Validation Cache: {stats['cache_sizes']['validation_cache']} entries
- Total: {stats['total_cached_entries']} entries

Cache Files:
- Total Size: {stats['total_cache_size_bytes']} bytes
- TTL: {stats['ttl_days']} days
- Last Cleanup: {stats['last_cleanup'] or 'Never'}
"""
        return info.strip()
