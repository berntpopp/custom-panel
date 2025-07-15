"""
Unit tests for gnomAD cache functionality.
"""

from __future__ import annotations

import tempfile
import time
from pathlib import Path
from unittest.mock import patch

import pytest

from custom_panel.core.gnomad_cache import GnomADCache


class TestGnomADCache:
    """Test the gnomAD cache functionality."""

    @pytest.fixture
    def temp_cache_dir(self):
        """Create a temporary cache directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            yield temp_dir

    @pytest.fixture
    def cache(self, temp_cache_dir):
        """Create a test cache instance."""
        return GnomADCache(temp_cache_dir, ttl_days=1)

    def test_cache_initialization(self, cache, temp_cache_dir):
        """Test cache initialization."""
        assert cache.cache_dir == Path(temp_cache_dir)
        assert cache.ttl_seconds == 24 * 60 * 60  # 1 day in seconds
        assert isinstance(cache.rsid_cache, dict)
        assert isinstance(cache.liftover_cache, dict)
        assert isinstance(cache.validation_cache, dict)

    def test_cache_key_generation(self, cache):
        """Test cache key generation methods."""
        # Test rsID cache key
        rsid_key = cache.get_rsid_cache_key("1", 12345, "A", "T", "GRCh38")
        assert rsid_key == "rsid:GRCh38:1:12345:A:T"

        # Test liftover cache key
        liftover_key = cache.get_liftover_cache_key("1-12345-A-T", "GRCh38", "GRCh37")
        assert liftover_key == "liftover:GRCh38:GRCh37:1-12345-A-T"

    def test_set_and_get_rsid_cache(self, cache):
        """Test setting and getting rsID cache entries."""
        key = "rsid:GRCh38:1:12345:A:T"
        value = "rs123456"

        # Set value
        cache.set(key, value)

        # Get value
        result = cache.get(key)
        assert result == value

    def test_set_and_get_liftover_cache(self, cache):
        """Test setting and getting liftover cache entries."""
        key = "liftover:GRCh38:GRCh37:1-12345-A-T"
        value = ("1", "1-54321-A-T", 54321, "A", "T")

        # Set value
        cache.set(key, value)

        # Get value
        result = cache.get(key)
        assert result == value

    def test_set_and_get_validation_cache(self, cache):
        """Test setting and getting validation cache entries."""
        key = "validate:GRCh38:1-12345-A-T"
        value = True

        # Set value
        cache.set(key, value)

        # Get value
        result = cache.get(key)
        assert result == value

    def test_get_unknown_key_prefix(self, cache):
        """Test getting value with unknown key prefix."""
        key = "unknown:some:key"
        result = cache.get(key)
        assert result is None

    def test_cache_expiration(self, cache):
        """Test cache entry expiration."""
        key = "rsid:GRCh38:1:12345:A:T"
        value = "rs123456"

        # Set value
        cache.set(key, value)

        # Mock time to simulate expiration
        with patch(
            "time.time", return_value=time.time() + 25 * 60 * 60
        ):  # 25 hours later
            result = cache.get(key)
            assert result is None

    def test_cache_persistence(self, temp_cache_dir):
        """Test cache persistence across instances."""
        # Create first cache instance and set value
        cache1 = GnomADCache(temp_cache_dir, ttl_days=1)
        key = "rsid:GRCh38:1:12345:A:T"
        value = "rs123456"
        cache1.set(key, value)

        # Create second cache instance and check if value persists
        cache2 = GnomADCache(temp_cache_dir, ttl_days=1)
        result = cache2.get(key)
        assert result == value

    def test_cleanup_expired_entries(self, cache):
        """Test cleanup of expired cache entries."""
        current_time = time.time()

        # Manually set entries with different timestamps
        cache.rsid_cache["rsid:GRCh38:1:12345:A:T"] = (
            "rs123456",
            current_time - 25 * 60 * 60,
        )  # Expired
        cache.rsid_cache["rsid:GRCh38:1:67890:G:C"] = (
            "rs789012",
            current_time,
        )  # Not expired

        expired_count = cache.cleanup_expired()
        assert expired_count == 1

        # Check that expired entry is gone
        assert "rsid:GRCh38:1:12345:A:T" not in cache.rsid_cache
        assert "rsid:GRCh38:1:67890:G:C" in cache.rsid_cache

    def test_clear_cache(self, cache):
        """Test clearing all cache entries."""
        # Set entries in different caches
        cache.set("rsid:GRCh38:1:12345:A:T", "rs123456")
        cache.set(
            "liftover:GRCh38:GRCh37:1-12345-A-T", ("1", "1-54321-A-T", 54321, "A", "T")
        )
        cache.set("validate:GRCh38:1-12345-A-T", True)

        # Clear cache
        cleared_count = cache.clear()
        assert cleared_count == 3

        # Check that all entries are gone
        assert cache.get("rsid:GRCh38:1:12345:A:T") is None
        assert cache.get("liftover:GRCh38:GRCh37:1-12345-A-T") is None
        assert cache.get("validate:GRCh38:1-12345-A-T") is None

    def test_cache_stats(self, cache):
        """Test cache statistics."""
        # Set some entries
        cache.set("rsid:GRCh38:1:12345:A:T", "rs123456")
        cache.set(
            "liftover:GRCh38:GRCh37:1-12345-A-T", ("1", "1-54321-A-T", 54321, "A", "T")
        )

        # Get entries (hits)
        cache.get("rsid:GRCh38:1:12345:A:T")
        cache.get("liftover:GRCh38:GRCh37:1-12345-A-T")

        # Get non-existent entry (miss)
        cache.get("rsid:GRCh38:1:99999:X:Y")

        stats = cache.get_stats()

        assert stats["hits"] == 2
        assert stats["misses"] == 1
        assert stats["sets"] == 2
        assert stats["total_requests"] == 3
        assert stats["hit_rate"] == 66.67
        assert stats["total_cached_entries"] == 2

    def test_cache_info_string(self, cache):
        """Test cache information string generation."""
        # Set some entries
        cache.set("rsid:GRCh38:1:12345:A:T", "rs123456")
        cache.get("rsid:GRCh38:1:12345:A:T")

        info = cache.get_cache_info()

        assert "gnomAD Cache Information" in info
        assert "Total Requests:" in info
        assert "Cache Hits:" in info
        assert "Cache Sizes:" in info

    def test_load_corrupted_cache_file(self, temp_cache_dir):
        """Test handling of corrupted cache files."""
        # Create a corrupted cache file
        cache_file = Path(temp_cache_dir) / "rsid_cache.pkl"
        cache_file.write_text("corrupted data")

        # Initialize cache - should handle corrupted file gracefully
        cache = GnomADCache(temp_cache_dir, ttl_days=1)

        # Should have empty cache
        assert len(cache.rsid_cache) == 0

    def test_save_cache_error_handling(self, cache):
        """Test error handling during cache save."""
        # Mock file write error
        with patch("builtins.open", side_effect=OSError("Write failed")):
            # This should not raise an exception
            cache.set("rsid:GRCh38:1:12345:A:T", "rs123456")

    def test_load_stats_error_handling(self, temp_cache_dir):
        """Test error handling during stats loading."""
        # Create corrupted stats file
        stats_file = Path(temp_cache_dir) / "cache_stats.json"
        stats_file.write_text("corrupted json")

        # Initialize cache - should handle corrupted stats gracefully
        cache = GnomADCache(temp_cache_dir, ttl_days=1)

        # Should have default stats
        assert cache.stats["hits"] == 0
        assert cache.stats["misses"] == 0
