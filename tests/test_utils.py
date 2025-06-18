"""Tests for utility functions."""


import pytest

from custom_panel.core.utils import get_source_group_from_config, normalize_count


class TestNormalizeCount:
    """Test normalization utility functions."""

    def test_logistic_normalization_basic(self):
        """Test basic logistic normalization."""
        params = {"k": 0.35, "x0": 5}

        # Test known points
        result = normalize_count(5, "logistic", params)
        assert abs(result - 0.5) < 0.01  # At midpoint should be ~0.5

        result = normalize_count(0, "logistic", params)
        assert result < 0.2  # Low count should be low confidence

        result = normalize_count(10, "logistic", params)
        assert result > 0.8  # High count should be high confidence

    def test_logistic_normalization_steep_curve(self):
        """Test logistic normalization with steep curve."""
        params = {"k": 1.0, "x0": 3}

        result = normalize_count(3, "logistic", params)
        assert abs(result - 0.5) < 0.01

        result = normalize_count(1, "logistic", params)
        assert result < 0.3

        result = normalize_count(5, "logistic", params)
        assert result > 0.7

    def test_linear_normalization_basic(self):
        """Test basic linear normalization."""
        params = {"max_count": 3}

        assert normalize_count(0, "linear", params) == 0.0
        assert normalize_count(1, "linear", params) == 1 / 3
        assert normalize_count(2, "linear", params) == 2 / 3
        assert normalize_count(3, "linear", params) == 1.0
        assert normalize_count(5, "linear", params) == 1.0  # Capped at 1.0

    def test_linear_normalization_large_max(self):
        """Test linear normalization with large max."""
        params = {"max_count": 10}

        assert normalize_count(0, "linear", params) == 0.0
        assert normalize_count(5, "linear", params) == 0.5
        assert normalize_count(10, "linear", params) == 1.0
        assert normalize_count(15, "linear", params) == 1.0

    def test_edge_cases(self):
        """Test edge cases and boundary conditions."""
        # Very large counts with logistic
        params = {"k": 0.35, "x0": 5}
        result = normalize_count(1000, "logistic", params)
        assert 0.99 <= result <= 1.0

        # Very small k value
        params = {"k": 0.01, "x0": 5}
        result = normalize_count(10, "logistic", params)
        assert 0.5 < result < 0.6  # Should be close to midpoint

        # Large k value (steep curve)
        params = {"k": 5.0, "x0": 3}
        result_low = normalize_count(2, "logistic", params)
        result_high = normalize_count(4, "logistic", params)
        assert result_low < 0.1
        assert result_high > 0.9

    def test_error_conditions(self):
        """Test error conditions."""
        # Negative count
        with pytest.raises(ValueError, match="Count must be non-negative"):
            normalize_count(-1, "logistic", {"k": 1, "x0": 5})

        # Missing parameters for logistic
        with pytest.raises(ValueError, match="Logistic normalization requires"):
            normalize_count(5, "logistic", {"k": 1})  # Missing x0

        with pytest.raises(ValueError, match="Logistic normalization requires"):
            normalize_count(5, "logistic", {"x0": 5})  # Missing k

        # Missing parameters for linear
        with pytest.raises(ValueError, match="Linear normalization requires"):
            normalize_count(5, "linear", {})

        # Invalid max_count for linear
        with pytest.raises(ValueError, match="max_count must be positive"):
            normalize_count(5, "linear", {"max_count": 0})

        with pytest.raises(ValueError, match="max_count must be positive"):
            normalize_count(5, "linear", {"max_count": -1})

        # Unknown method
        with pytest.raises(ValueError, match="Unknown normalization method"):
            normalize_count(5, "unknown_method", {})

    def test_overflow_protection(self):
        """Test protection against mathematical overflow."""
        # Test extreme parameter values that could cause overflow
        params = {"k": 100, "x0": 5}

        # Should not raise overflow error
        result_low = normalize_count(0, "logistic", params)
        result_high = normalize_count(10, "logistic", params)

        assert 0 <= result_low <= 1
        assert 0 <= result_high <= 1


class TestGetSourceGroupFromConfig:
    """Test source group configuration utilities."""

    def test_source_group_mapping(self):
        """Test mapping sources to source groups."""
        config = {
            "data_sources": {
                "Commercial_Panels": {
                    "source_group": True,
                    "panels": [
                        {"name": "Myriad_myRisk"},
                        {"name": "Blueprint_Genetics"},
                    ],
                },
                "PanelApp": {"source_group": True, "panels": [{"name": "PanelApp_UK"}]},
                "ACMG_Incidental_Findings": {"source_group": False},
            }
        }

        # Test source group members
        assert (
            get_source_group_from_config("Myriad_myRisk", config) == "Commercial_Panels"
        )
        assert (
            get_source_group_from_config("Blueprint_Genetics", config)
            == "Commercial_Panels"
        )
        assert get_source_group_from_config("PanelApp_UK", config) == "PanelApp"

        # Test standalone source
        assert (
            get_source_group_from_config("ACMG_Incidental_Findings", config)
            == "ACMG_Incidental_Findings"
        )

        # Test unknown source
        assert get_source_group_from_config("Unknown_Source", config) is None

    def test_empty_config(self):
        """Test with empty or missing configuration."""
        assert get_source_group_from_config("Any_Source", {}) is None
        assert get_source_group_from_config("Any_Source", {"data_sources": {}}) is None

    def test_malformed_config(self):
        """Test with malformed configuration."""
        config = {
            "data_sources": {
                "Bad_Group": "not_a_dict",
                "Good_Group": {
                    "source_group": True,
                    "panels": [{"name": "Test_Source"}],
                },
            }
        }

        assert get_source_group_from_config("Test_Source", config) == "Good_Group"
        assert get_source_group_from_config("Bad_Group", config) is None
