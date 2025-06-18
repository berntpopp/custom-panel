"""
Tests for CLI integration with commercial panels.

This module tests that the CLI properly integrates the commercial panels source.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import patch

from custom_panel.cli import fetch_all_sources


class TestCLIIntegration:
    """Test CLI integration with commercial panels."""

    def test_fetch_all_sources_includes_commercial_panels(self):
        """Test that fetch_all_sources includes commercial panels when enabled."""

        # Create a temporary JSON file for testing
        with tempfile.TemporaryDirectory() as temp_dir:
            json_path = Path(temp_dir) / "test_panel.json"

            test_data = {
                "panel_name": "Test Panel",
                "source_url": "http://test.com",
                "retrieval_date": "2024-01-01",
                "genes": ["BRCA1", "BRCA2"],
            }

            with open(json_path, "w") as f:
                json.dump(test_data, f)

            # Create config with commercial panels enabled
            config = {
                "data_sources": {
                    "PanelApp": {"enabled": False},
                    "Inhouse_Panels": {"enabled": False},
                    "ACMG_Incidental_Findings": {"enabled": False},
                    "Manual_Curation": {"enabled": False},
                    "HPO_Neoplasm": {"enabled": False},
                    "COSMIC_Germline": {"enabled": False},
                    "Commercial_Panels": {
                        "enabled": True,
                        "source_group": True,
                        "panels": [
                            {
                                "name": "Test Panel",
                                "file_path": str(json_path),
                                "evidence_score": 0.8,
                            }
                        ],
                    },
                }
            }

            # Mock the console to avoid output during testing
            with patch("custom_panel.cli.console"):
                dataframes = fetch_all_sources(config)

            # Should have one DataFrame from commercial panels
            assert len(dataframes) == 1
            df = dataframes[0]

            # Should have the expected genes
            assert len(df) == 2
            assert "BRCA1" in df["approved_symbol"].values
            assert "BRCA2" in df["approved_symbol"].values

    def test_commercial_panels_disabled(self):
        """Test that commercial panels are skipped when disabled."""
        config = {
            "data_sources": {
                "PanelApp": {"enabled": False},
                "Inhouse_Panels": {"enabled": False},
                "ACMG_Incidental_Findings": {"enabled": False},
                "Manual_Curation": {"enabled": False},
                "HPO_Neoplasm": {"enabled": False},
                "COSMIC_Germline": {"enabled": False},
                "Commercial_Panels": {"enabled": False},
            }
        }

        # Mock the console to avoid output during testing
        with patch("custom_panel.cli.console"):
            dataframes = fetch_all_sources(config)

        # Should have no DataFrames since all sources are disabled
        assert len(dataframes) == 0

    def test_fetch_command_includes_commercial(self):
        """Test that the fetch command handles commercial source."""
        from custom_panel.cli import fetch_all_sources

        # Test that the fetch_all_sources function has commercial_panels in source_functions
        # This verifies the integration was done correctly
        config = {
            "data_sources": {"Commercial_Panels": {"enabled": True, "panels": []}}
        }

        # This should not error out - if commercial_panels wasn't integrated,
        # it would raise a KeyError
        with patch("custom_panel.cli.console"):
            dataframes = fetch_all_sources(config)

        # Should work without error (empty result since no panels configured)
        assert isinstance(dataframes, list)
