"""
Unit tests for gnomAD client functionality.
"""

from __future__ import annotations

from unittest.mock import Mock, patch

import pytest
import requests

from custom_panel.core.gnomad_client import GnomADAPIError, GnomADClient, RateLimiter


class TestRateLimiter:
    """Test the rate limiter functionality."""

    def test_rate_limiter_initialization(self):
        """Test rate limiter initialization."""
        limiter = RateLimiter(5.0)
        assert limiter.requests_per_second == 5.0
        assert limiter.min_interval == 0.2

    @patch("time.time")
    @patch("time.sleep")
    def test_rate_limiter_enforces_limit(self, mock_sleep, mock_time):
        """Test that rate limiter enforces the rate limit."""
        mock_time.side_effect = [
            0.1,
            0.1,
        ]  # Current time is 0.1, last request was at 0.0

        limiter = RateLimiter(5.0)  # 5 requests per second = 0.2s minimum interval
        limiter.last_request_time = 0.0

        limiter.wait_if_needed()

        mock_sleep.assert_called_once_with(
            0.1
        )  # Should sleep for 0.1s to reach 0.2s interval


class TestGnomADClient:
    """Test the gnomAD client functionality."""

    @pytest.fixture
    def client(self):
        """Create a test gnomAD client."""
        with patch("custom_panel.core.gnomad_client.GnomADCache"):
            return GnomADClient(rate_limit=10.0, cache_ttl_days=1)

    def test_client_initialization(self, client):
        """Test client initialization."""
        assert client.base_url == "https://gnomad.broadinstitute.org/api"
        assert client.timeout == 30
        assert client.retry_attempts == 3

    @patch("requests.Session.post")
    def test_make_graphql_request_success(self, mock_post, client):
        """Test successful GraphQL request."""
        # Mock successful response
        mock_response = Mock()
        mock_response.json.return_value = {"data": {"variant": {"rsid": "rs123456"}}}
        mock_response.raise_for_status.return_value = None
        mock_post.return_value = mock_response

        # Test request
        query = "query { variant { rsid } }"
        variables = {"variantId": "1-12345-A-T"}

        result = client._make_graphql_request(query, variables)

        assert result == {"variant": {"rsid": "rs123456"}}
        mock_post.assert_called_once()

    @patch("requests.Session.post")
    def test_make_graphql_request_graphql_error(self, mock_post, client):
        """Test GraphQL error handling."""
        # Mock GraphQL error response
        mock_response = Mock()
        mock_response.json.return_value = {
            "errors": [{"message": "Variable not found"}]
        }
        mock_response.raise_for_status.return_value = None
        mock_post.return_value = mock_response

        query = "query { variant { rsid } }"
        variables = {"variantId": "1-12345-A-T"}

        with pytest.raises(GnomADAPIError, match="GraphQL errors"):
            client._make_graphql_request(query, variables)

    @patch("requests.Session.post")
    def test_make_graphql_request_network_error(self, mock_post, client):
        """Test network error handling."""
        # Mock network error
        mock_post.side_effect = requests.exceptions.RequestException("Network error")

        query = "query { variant { rsid } }"
        variables = {"variantId": "1-12345-A-T"}

        with pytest.raises(GnomADAPIError, match="Request failed"):
            client._make_graphql_request(query, variables)

    def test_resolve_rsid_cache_hit(self, client):
        """Test rsID resolution with cache hit."""
        # Mock cache hit
        client.cache.get.return_value = "rs123456"

        result = client.resolve_rsid("1", 12345, "A", "T", "GRCh38")

        assert result == "rs123456"
        client.cache.get.assert_called_once()

    @patch.object(GnomADClient, "_make_graphql_request")
    def test_resolve_rsid_success(self, mock_request, client):
        """Test successful rsID resolution."""
        # Mock cache miss
        client.cache.get.return_value = None

        # Mock successful API response
        mock_request.return_value = {
            "variant": {"variantId": "1-12345-A-T", "rsid": "rs123456"}
        }

        result = client.resolve_rsid("1", 12345, "A", "T", "GRCh38")

        assert result == "rs123456"
        client.cache.set.assert_called_once()

    @patch.object(GnomADClient, "_make_graphql_request")
    def test_resolve_rsid_not_found(self, mock_request, client):
        """Test rsID resolution when variant not found."""
        # Mock cache miss
        client.cache.get.return_value = None

        # Mock API response with no variant
        mock_request.return_value = {"variant": None}

        result = client.resolve_rsid("1", 12345, "A", "T", "GRCh38")

        assert result is None

    def test_liftover_coordinates_same_build(self, client):
        """Test liftover with same source and target build."""
        result = client.liftover_coordinates("1-12345-A-T", "GRCh38", "GRCh38")

        assert result == ("1", "1-12345-A-T", 12345, "A", "T")

    @patch.object(GnomADClient, "_make_graphql_request")
    def test_liftover_coordinates_success(self, mock_request, client):
        """Test successful coordinate liftover."""
        # Mock cache miss
        client.cache.get.return_value = None

        # Mock successful liftover response
        mock_request.return_value = {
            "liftover": {
                "source": {"variant_id": "1-12345-A-T", "reference_genome": "GRCh38"},
                "liftover": {"variant_id": "1-54321-A-T", "reference_genome": "GRCh37"},
            }
        }

        result = client.liftover_coordinates("1-12345-A-T", "GRCh38", "GRCh37")

        assert result == ("1", "1-54321-A-T", 54321, "A", "T")

    @patch.object(GnomADClient, "_make_graphql_request")
    def test_validate_variant_success(self, mock_request, client):
        """Test successful variant validation."""
        # Mock cache miss
        client.cache.get.return_value = None

        # Mock successful validation response
        mock_request.return_value = {
            "variant": {
                "variantId": "1-12345-A-T",
                "chrom": "1",
                "pos": 12345,
                "ref": "A",
                "alt": "T",
            }
        }

        result = client.validate_variant("1", 12345, "A", "T", "GRCh38")

        assert result is True

    @patch.object(GnomADClient, "_make_graphql_request")
    def test_validate_variant_not_found(self, mock_request, client):
        """Test variant validation when variant not found."""
        # Mock cache miss
        client.cache.get.return_value = None

        # Mock validation response with no variant
        mock_request.return_value = {"variant": None}

        result = client.validate_variant("1", 12345, "A", "T", "GRCh38")

        assert result is False

    def test_get_cache_stats(self, client):
        """Test getting cache statistics."""
        client.cache.get_stats.return_value = {"hits": 10, "misses": 5}

        stats = client.get_cache_stats()

        assert stats == {"hits": 10, "misses": 5}

    def test_clear_cache(self, client):
        """Test clearing cache."""
        client.cache.clear.return_value = 15

        result = client.clear_cache()

        assert result == 15

    @patch("requests.Session.post")
    def test_graphql_request_429_retry_success(self, mock_post, client):
        """Test successful retry after 429 rate limit error."""
        # First response: 429 rate limit
        mock_429_response = Mock()
        mock_429_response.status_code = 429
        mock_429_response.reason = "Too Many Requests"

        # Second response: success
        mock_success_response = Mock()
        mock_success_response.status_code = 200
        mock_success_response.json.return_value = {
            "data": {"variant": {"rsid": "rs123456"}}
        }
        mock_success_response.raise_for_status.return_value = None

        mock_post.side_effect = [mock_429_response, mock_success_response]

        result = client._make_graphql_request("query test", {})

        assert result == {"variant": {"rsid": "rs123456"}}
        assert mock_post.call_count == 2

    @patch("requests.Session.post")
    def test_graphql_request_429_max_retries_exceeded(self, mock_post, client):
        """Test failure after max 429 retries exceeded."""
        mock_response = Mock()
        mock_response.status_code = 429
        mock_response.reason = "Too Many Requests"
        mock_post.return_value = mock_response

        with pytest.raises(GnomADAPIError) as exc_info:
            client._make_graphql_request("query test", {}, max_429_retries=2)

        assert "too many 429 error responses after 2 retries" in str(exc_info.value)
        assert mock_post.call_count == 3  # Initial + 2 retries

    def test_rate_limiter_exponential_backoff(self, client):
        """Test rate limiter exponential backoff for 429 errors."""
        # Test initial state
        assert client.rate_limiter.consecutive_429_errors == 0

        # Record 429 errors and verify backoff increases
        client.rate_limiter.record_429_error()
        assert client.rate_limiter.consecutive_429_errors == 1

        client.rate_limiter.record_429_error()
        assert client.rate_limiter.consecutive_429_errors == 2

        # Reset and verify
        client.rate_limiter.reset_429_errors()
        assert client.rate_limiter.consecutive_429_errors == 0

    @patch("time.sleep")
    @patch("time.time")
    def test_rate_limiter_exponential_backoff_timing(
        self, mock_time, mock_sleep, client
    ):
        """Test rate limiter calculates exponential backoff delays correctly."""
        mock_time.return_value = 0.0

        limiter = client.rate_limiter
        limiter.last_request_time = 0.0

        # First 429 error - should use 1s base delay
        limiter.record_429_error()
        limiter.wait_if_needed()

        # Should sleep for at least the base 429 delay
        assert mock_sleep.call_count >= 1

        # Reset for next test
        mock_sleep.reset_mock()

        # Second 429 error - should use 2s (base * 2^1)
        limiter.record_429_error()
        limiter.wait_if_needed()

        # Should sleep for exponentially increased time
        assert mock_sleep.call_count >= 1

    @patch("requests.Session.post")
    def test_graphql_request_connection_pool_429_retry(self, mock_post, client):
        """Test retry on connection pool 429 error."""
        # First request: connection pool error with 429 mention
        mock_post.side_effect = [
            requests.exceptions.RequestException("too many 429 error responses"),
            Mock(
                status_code=200,
                json=lambda: {"data": {"variant": {"rsid": "rs123456"}}},
                raise_for_status=lambda: None,
            ),
        ]

        result = client._make_graphql_request("query test", {}, max_429_retries=2)

        assert result == {"variant": {"rsid": "rs123456"}}
        assert mock_post.call_count == 2
