"""
Tests for CLI integration with commercial panels.

This module tests that the CLI properly integrates the commercial panels source.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import patch

from custom_panel.engine.pipeline import Pipeline


class TestCLIIntegration:
    """Test CLI integration with commercial panels."""

    def test_pipeline_includes_commercial_panels(self):
        """Test that pipeline includes commercial panels when enabled."""

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
                "general": {"output_dir": str(temp_dir)},
                "data_sources": {
                    "PanelApp": {"enabled": False},
                    "Inhouse_Panels": {"enabled": False},
                    "ACMG_Incidental_Findings": {"enabled": False},
                    "Manual_Curation": {"enabled": False},
                    "HPO_Neoplasm": {"enabled": False},
                    "COSMIC_Germline": {"enabled": False},
                    "ClinGen": {"enabled": False},
                    "TheGenCC": {"enabled": False},
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
                },
            }

            # Create pipeline and test source fetching
            pipeline = Pipeline(config)

            # Mock the console to avoid output during testing
            with patch("custom_panel.engine.pipeline.console"):
                dataframes = pipeline._fetch_all_sources()

            # Should have one DataFrame from commercial panels
            assert len(dataframes) == 1
            df = dataframes[0]

            # Should have the expected genes
            assert len(df) == 2
            assert "BRCA1" in df["approved_symbol"].values
            assert "BRCA2" in df["approved_symbol"].values

    def test_commercial_panels_disabled(self):
        """Test that commercial panels are skipped when disabled."""
        with tempfile.TemporaryDirectory() as temp_dir:
            config = {
                "general": {"output_dir": str(temp_dir)},
                "data_sources": {
                    "PanelApp": {"enabled": False},
                    "Inhouse_Panels": {"enabled": False},
                    "ACMG_Incidental_Findings": {"enabled": False},
                    "Manual_Curation": {"enabled": False},
                    "HPO_Neoplasm": {"enabled": False},
                    "COSMIC_Germline": {"enabled": False},
                    "ClinGen": {"enabled": False},
                    "TheGenCC": {"enabled": False},
                    "Commercial_Panels": {"enabled": False},
                },
            }

            # Create pipeline and test source fetching
            pipeline = Pipeline(config)

            # Mock the console to avoid output during testing
            with patch("custom_panel.engine.pipeline.console"):
                dataframes = pipeline._fetch_all_sources()

            # Should have no DataFrames since all sources are disabled
            assert len(dataframes) == 0

    def test_pipeline_includes_commercial_source(self):
        """Test that the pipeline handles commercial source correctly."""

        with tempfile.TemporaryDirectory() as temp_dir:
            # Test that the pipeline has commercial_panels in source_functions
            # This verifies the integration was done correctly
            config = {
                "general": {"output_dir": str(temp_dir)},
                "data_sources": {"Commercial_Panels": {"enabled": True, "panels": []}},
            }

            # Create pipeline and test source fetching
            pipeline = Pipeline(config)

            # This should not error out - if commercial_panels wasn't integrated,
            # it would raise a KeyError
            with patch("custom_panel.engine.pipeline.console"):
                dataframes = pipeline._fetch_all_sources()

            # Should work without error (empty result since no panels configured)
            assert isinstance(dataframes, list)
