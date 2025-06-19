"""
Tests for the scrapers framework.

This module tests the individual parser classes and the master runner script.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from scrapers.parsers.parse_blueprint import BlueprintParser
from scrapers.parsers.parse_fulgent import FulgentParser
from scrapers.parsers.parse_myriad import MyriadParser
from scrapers.run_scrapers import create_output_json, get_parser_class, run_scraper


class TestParsers:
    """Test individual parser classes."""

    def test_myriad_parser(self):
        """Test MyriadParser with mock HTML content."""
        # Load fixture data
        fixture_path = (
            Path(__file__).parent / "data" / "scraper_fixtures" / "myriad_myrisk.html"
        )
        with open(fixture_path) as f:
            mock_html = f.read()

        # Mock the session.get call
        mock_response = MagicMock()
        mock_response.content = mock_html.encode("utf-8")
        mock_response.raise_for_status.return_value = None

        mock_session = MagicMock()
        mock_session.get.return_value = mock_response

        with patch("requests.Session", return_value=mock_session):
            parser = MyriadParser("http://test.com")
            genes = parser.parse()

        # Expected genes from the fixture
        expected_genes = ["BRCA1", "BRCA2", "TP53", "MLH1", "MSH2"]

        assert len(genes) == len(expected_genes)
        for gene in expected_genes:
            assert gene in genes

    def test_blueprint_parser(self):
        """Test BlueprintParser with mock HTML content."""
        # Load fixture data
        fixture_path = (
            Path(__file__).parent
            / "data"
            / "scraper_fixtures"
            / "blueprint_genetics.html"
        )
        with open(fixture_path) as f:
            mock_html = f.read()

        # Mock the requests.get call
        mock_response = MagicMock()
        mock_response.content = mock_html.encode("utf-8")
        mock_response.raise_for_status.return_value = None

        with patch("requests.get", return_value=mock_response):
            parser = BlueprintParser("http://test.com")
            genes = parser.parse()

        # Expected genes from the fixture
        expected_genes = ["APC", "ATM", "CHEK2", "PALB2"]

        assert len(genes) == len(expected_genes)
        for gene in expected_genes:
            assert gene in genes

    def test_fulgent_parser(self):
        """Test FulgentParser with mock HTML content."""
        # Load fixture data
        fixture_path = (
            Path(__file__).parent
            / "data"
            / "scraper_fixtures"
            / "fulgent_genetics.html"
        )
        with open(fixture_path) as f:
            mock_html = f.read()

        # Mock the requests.get call
        mock_response = MagicMock()
        mock_response.content = mock_html.encode("utf-8")
        mock_response.raise_for_status.return_value = None

        with patch("requests.get", return_value=mock_response):
            parser = FulgentParser("http://test.com")
            genes = parser.parse()

        # Expected genes from the meta description
        expected_genes = [
            "AIP",
            "BRCA1",
            "BRCA2",
            "TP53",
            "PTEN",
            "STK11",
            "CDH1",
            "PALB2",
            "ATM",
            "CHEK2",
            "NBN",
        ]

        assert len(genes) == len(expected_genes)
        for gene in expected_genes:
            assert gene in genes

    def test_gene_symbol_cleaning(self):
        """Test gene symbol cleaning and validation."""
        parser = MyriadParser("http://test.com")

        # Test cleaning
        assert parser.clean_gene_symbol("BRCA1*") == "BRCA1"
        assert parser.clean_gene_symbol("BRCA2 (some info)") == "BRCA2"
        assert parser.clean_gene_symbol("  TP53  ") == "TP53"
        assert parser.clean_gene_symbol("gene1,gene2") == "GENE1"

        # Test validation
        assert parser.validate_gene_symbol("BRCA1") is True
        assert parser.validate_gene_symbol("TP53") is True
        assert parser.validate_gene_symbol("ATM") is True
        assert parser.validate_gene_symbol("A") is True  # Minimum length
        assert (
            parser.validate_gene_symbol("ABCDEFGHIJ1234567890") is True
        )  # Maximum length

        # Test invalid symbols
        assert parser.validate_gene_symbol("") is False
        assert parser.validate_gene_symbol("123") is False  # No letters
        assert parser.validate_gene_symbol("GENE_WITH_SPECIAL@CHARS") is False
        assert parser.validate_gene_symbol("GENE") is False  # Skip term
        assert parser.validate_gene_symbol("DNA") is False  # Skip term


class TestRunnerScript:
    """Test the master runner script functionality."""

    def test_get_parser_class(self):
        """Test dynamic parser class loading."""
        # Test successful import
        parser_cls = get_parser_class("parse_myriad", "MyriadParser")
        assert parser_cls == MyriadParser

        # Test failed import
        with pytest.raises(ImportError):
            get_parser_class("nonexistent_module", "NonexistentParser")

    def test_create_output_json(self):
        """Test JSON output creation."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test_output.json"

            create_output_json(
                panel_name="Test Panel",
                source_url="http://test.com",
                genes=["BRCA1", "BRCA2", "TP53"],
                output_path=output_path,
            )

            # Verify file was created
            assert output_path.exists()

            # Verify content
            with open(output_path) as f:
                data = json.load(f)

            assert data["panel_name"] == "Test Panel"
            assert data["source_url"] == "http://test.com"
            assert "retrieval_date" in data
            assert data["genes"] == ["BRCA1", "BRCA2", "TP53"]

    def test_run_scraper_single_url(self):
        """Test running a single scraper with one URL."""
        scraper_config = {
            "url": "http://test.com",
            "parser_module": "parse_myriad",
            "parser_class": "MyriadParser",
            "output_path": "test_output.json",
        }

        # Mock the parser
        mock_parser_instance = MagicMock()
        mock_parser_instance.parse.return_value = ["BRCA1", "BRCA2"]

        mock_parser_class = MagicMock()
        mock_parser_class.return_value = mock_parser_instance

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test_output.json"
            scraper_config["output_path"] = str(output_path)

            with patch(
                "scrapers.run_scrapers.get_parser_class", return_value=mock_parser_class
            ):
                run_scraper("test_scraper", scraper_config)

            # Verify output file was created
            assert output_path.exists()

            # Verify content
            with open(output_path) as f:
                data = json.load(f)

            assert data["panel_name"] == "test_scraper"
            assert data["genes"] == ["BRCA1", "BRCA2"]

    def test_run_scraper_multiple_urls(self):
        """Test running a scraper with multiple URLs (like Blueprint)."""
        scraper_config = {
            "url": ["http://test1.com", "http://test2.com"],
            "parser_module": "parse_blueprint",
            "parser_class": "BlueprintParser",
            "output_path": "test_output.json",
        }

        # Mock the parser to return different genes for each URL
        mock_parser_instance = MagicMock()
        mock_parser_instance.parse.side_effect = [
            ["BRCA1", "BRCA2"],  # First URL
            ["TP53", "ATM"],  # Second URL
        ]

        mock_parser_class = MagicMock()
        mock_parser_class.return_value = mock_parser_instance

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test_output.json"
            scraper_config["output_path"] = str(output_path)

            with patch(
                "scrapers.run_scrapers.get_parser_class", return_value=mock_parser_class
            ):
                run_scraper("test_scraper", scraper_config)

            # Verify output file was created
            assert output_path.exists()

            # Verify content - should have genes from both URLs
            with open(output_path) as f:
                data = json.load(f)

            assert data["panel_name"] == "test_scraper"
            assert "BRCA1" in data["genes"]
            assert "BRCA2" in data["genes"]
            assert "TP53" in data["genes"]
            assert "ATM" in data["genes"]


class TestCommercialPanelsIntegration:
    """Test integration with the commercial panels source module."""

    def test_commercial_panels_source_module(self):
        """Test that the commercial panels source module works correctly."""
        from custom_panel.sources.g03_commercial_panels import (
            fetch_commercial_panels_data,
        )

        # Create a mock config
        config = {
            "data_sources": {
                "commercial_panels": {
                    "enabled": True,
                    "panels": [
                        {
                            "name": "Test Panel",
                            "file_path": "nonexistent.json",
                            "evidence_score": 0.8,
                        }
                    ],
                }
            }
        }

        # Should return empty DataFrame for nonexistent file
        df = fetch_commercial_panels_data(config)
        assert df.empty

    def test_commercial_panels_with_valid_json(self):
        """Test commercial panels with a valid JSON file."""
        from custom_panel.sources.g03_commercial_panels import (
            fetch_commercial_panels_data,
        )

        # Create a temporary JSON file
        with tempfile.TemporaryDirectory() as temp_dir:
            json_path = Path(temp_dir) / "test_panel.json"

            test_data = {
                "panel_name": "Test Panel",
                "source_url": "http://test.com",
                "retrieval_date": "2024-01-01",
                "genes": ["BRCA1", "BRCA2", "TP53"],
            }

            with open(json_path, "w") as f:
                json.dump(test_data, f)

            config = {
                "data_sources": {
                    "Commercial_Panels": {
                        "enabled": True,
                        "panels": [
                            {
                                "name": "Test Panel",
                                "file_path": str(json_path),
                                "evidence_score": 0.8,
                            }
                        ],
                    }
                }
            }

            df = fetch_commercial_panels_data(config)

            # Should have 3 genes
            assert len(df) == 3
            assert "BRCA1" in df["approved_symbol"].values
            assert "BRCA2" in df["approved_symbol"].values
            assert "TP53" in df["approved_symbol"].values

            # Check evidence scores
            assert all(score == 0.8 for score in df["source_evidence_score"].values)

            # Check source details
            assert all(
                "URL:http://test.com" in detail
                for detail in df["source_details"].values
            )
            assert all(
                "Date:2024-01-01" in detail for detail in df["source_details"].values
            )
