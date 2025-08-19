"""
Tests for data source extractors.
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd

from custom_panel.sources.a_incidental_findings import (
    DEFAULT_ACMG_GENES,
    _scrape_acmg_genes_from_ncbi,
    fetch_acmg_incidental_data,
)
from custom_panel.sources.b_manual_curation import (
    fetch_manual_curation_data,
    get_manual_curation_summary,
    process_manual_list,
    validate_manual_curation_config,
)
from custom_panel.sources.g00_inhouse_panels import (
    extract_genes_from_column,
    fetch_inhouse_panels_data,
    process_inhouse_panel,
    validate_inhouse_panel_config,
)
from custom_panel.sources.g01_panelapp import PanelAppClient, fetch_panelapp_data
from custom_panel.sources.g04_cosmic_germline import (
    _calculate_cosmic_score,
    _is_cache_valid,
    fetch_cosmic_germline_data,
    get_cosmic_summary,
    validate_cosmic_config,
)


class TestPanelAppClient:
    """Test PanelApp client functionality."""

    def test_init(self):
        """Test client initialization."""
        client = PanelAppClient(timeout=10, max_retries=2)
        assert client.timeout == 10
        assert client.max_retries == 2

    @patch("custom_panel.sources.g01_panelapp.requests.Session.get")
    def test_get_panel_list(self, mock_get):
        """Test getting panel list."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "results": [
                {"id": 1, "name": "Test Panel 1"},
                {"id": 2, "name": "Test Panel 2"},
            ],
            "next": None,
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = PanelAppClient()
        panels = client.get_panel_list("https://panelapp.example.com/api/v1/panels")

        assert len(panels) == 2
        assert panels[0]["id"] == 1
        assert panels[0]["name"] == "Test Panel 1"

    @patch("custom_panel.sources.g01_panelapp.requests.Session.get")
    def test_get_panel_genes(self, mock_get):
        """Test getting genes for a panel."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "results": [
                {
                    "gene_data": {"gene_symbol": "BRCA1"},
                    "confidence_level": "Green",
                    "mode_of_inheritance": "Autosomal dominant",
                    "phenotypes": [{"phenotype": "Breast cancer"}],
                    "publications": [{"pmid": 12345}],
                },
                {
                    "gene_data": {"gene_symbol": "TP53"},
                    "confidence_level": "Amber",
                    "mode_of_inheritance": "Autosomal dominant",
                    "phenotypes": [],
                    "publications": [],
                },
            ],
            "next": None,
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = PanelAppClient()
        genes = client.get_panel_genes("https://panelapp.example.com/api/v1/panels", 1)

        assert len(genes) == 2
        assert genes[0]["gene_data"]["gene_symbol"] == "BRCA1"
        assert genes[0]["confidence_level"] == "Green"


class TestFetchPanelAppData:
    """Test PanelApp data fetching."""

    @patch("custom_panel.sources.g01_panelapp.PanelAppClient")
    def test_fetch_panelapp_data(self, mock_client_class):
        """Test fetching PanelApp data with configuration."""
        # Mock client
        mock_client = Mock()
        mock_client_class.return_value = mock_client

        # Mock panel genes response
        mock_client.get_panel_genes.return_value = [
            {
                "gene_data": {"gene_symbol": "BRCA1"},
                "confidence_level": "Green",
                "mode_of_inheritance": "Autosomal dominant",
                "phenotypes": [{"phenotype": "Breast cancer"}],
                "publications": [{"pmid": 12345}],
            }
        ]

        config = {
            "data_sources": {
                "PanelApp": {
                    "enabled": True,
                    "panels": [
                        {
                            "name": "TestApp",
                            "base_url": "https://test.panelapp.com/api/v1/panels",
                            "panels": [{"id": 1, "name": "Test Panel"}],
                        }
                    ],
                    "evidence_scores": {"Green": 1.0, "Amber": 0.5, "Red": 0.1},
                }
            }
        }

        df = fetch_panelapp_data(config)

        assert len(df) == 1
        assert df.iloc[0]["approved_symbol"] == "BRCA1"
        assert df.iloc[0]["source_evidence_score"] == 1.0
        assert "TestApp:Test Panel" in df.iloc[0]["source_name"]

    def test_fetch_panelapp_data_disabled(self):
        """Test when PanelApp is disabled."""
        config = {"data_sources": {"PanelApp": {"enabled": False}}}

        df = fetch_panelapp_data(config)
        assert df.empty


class TestInhousePanels:
    """Test in-house panels functionality."""

    def test_extract_genes_from_column(self):
        """Test extracting genes from DataFrame column."""
        df = pd.DataFrame(
            {
                "Gene": ["BRCA1", "TP53,EGFR", "KRAS"],
                "Other": ["data1", "data2", "data3"],
            }
        )

        genes = extract_genes_from_column(df, "Gene")

        assert len(genes) == 4  # BRCA1, TP53, EGFR, KRAS
        assert "BRCA1" in genes
        assert "TP53" in genes
        assert "EGFR" in genes
        assert "KRAS" in genes

    def test_extract_genes_semicolon_separated(self):
        """Test extracting semicolon-separated genes."""
        df = pd.DataFrame({"Genes": ["BRCA1;TP53", "EGFR;KRAS;PIK3CA"]})

        genes = extract_genes_from_column(df, "Genes")

        assert len(genes) == 5
        assert all(
            gene in genes for gene in ["BRCA1", "TP53", "EGFR", "KRAS", "PIK3CA"]
        )

    def test_extract_genes_missing_column(self):
        """Test extracting genes from missing column."""
        df = pd.DataFrame({"Other": ["BRCA1", "TP53"]})

        genes = extract_genes_from_column(df, "Gene")
        assert genes == []

    def test_validate_inhouse_panel_config_valid(self):
        """Test validation of valid config."""
        config = {
            "name": "Test Panel",
            "file_path": __file__,  # This file exists
            "gene_column": "Gene",
            "evidence_score": 1.0,
        }

        errors = validate_inhouse_panel_config(config)
        assert len(errors) == 0

    def test_validate_inhouse_panel_config_invalid(self):
        """Test validation of invalid config."""
        config = {"file_path": "/nonexistent/file.xlsx", "evidence_score": -1.0}

        errors = validate_inhouse_panel_config(config)
        assert len(errors) > 0
        assert any("name is required" in error for error in errors)
        assert any("File not found" in error for error in errors)
        assert any("Evidence score must be" in error for error in errors)

    def test_process_inhouse_panel_csv(self):
        """Test processing CSV in-house panel."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("Gene Symbol,Description\n")
            f.write("BRCA1,Breast cancer gene\n")
            f.write("TP53,Tumor suppressor\n")
            f.flush()

            try:
                config = {
                    "name": "Test_CSV_Panel",
                    "file_path": f.name,
                    "gene_column": "Gene Symbol",
                    "evidence_score": 0.9,
                }

                df = process_inhouse_panel(config)

                assert df is not None
                assert len(df) == 2
                assert "BRCA1" in df["approved_symbol"].values
                assert "TP53" in df["approved_symbol"].values
                assert all(df["source_evidence_score"] == 0.9)
                assert all(df["source_name"] == "Test_CSV_Panel")
                # Test default category behavior
                assert "category" in df.columns
                assert all(df["category"] == "germline")  # Should default to germline
                assert all(
                    "Category:germline" in detail for detail in df["source_details"]
                )
            finally:
                Path(f.name).unlink()

    def test_process_inhouse_panel_text(self):
        """Test processing text file in-house panel."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\n")
            f.write("TP53\n")
            f.write("EGFR\n")
            f.flush()

            try:
                config = {
                    "name": "Test_Text_Panel",
                    "file_path": f.name,
                    "gene_column": "Gene",  # Will be used as column name
                    "evidence_score": 1.0,
                }

                df = process_inhouse_panel(config)

                assert df is not None
                assert len(df) == 3
                assert set(df["approved_symbol"]) == {"BRCA1", "TP53", "EGFR"}
            finally:
                Path(f.name).unlink()

    def test_fetch_inhouse_panels_data(self):
        """Test fetching in-house panels data."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\nTP53\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "Inhouse_Panels": {
                            "enabled": True,
                            "panels": [
                                {
                                    "name": "Test_Panel",
                                    "file_path": f.name,
                                    "gene_column": "Gene",
                                    "evidence_score": 1.0,
                                }
                            ],
                        }
                    }
                }

                df = fetch_inhouse_panels_data(config)

                assert len(df) == 2
                assert "BRCA1" in df["approved_symbol"].values
                assert "TP53" in df["approved_symbol"].values
            finally:
                Path(f.name).unlink()


class TestACMGIncidentalFindings:
    """Test ACMG incidental findings functionality."""

    def test_fetch_acmg_incidental_data_default(self):
        """Test fetching ACMG data with default gene list."""
        config = {
            "data_sources": {
                "acmg_incidental": {"enabled": True, "evidence_score": 1.5}
            }
        }

        df = fetch_acmg_incidental_data(config)

        assert len(df) == len(DEFAULT_ACMG_GENES)
        assert all(df["source_evidence_score"] == 1.5)
        assert all(df["source_name"] == "ACMG_Incidental_Findings")
        assert "BRCA1" in df["approved_symbol"].values
        assert "BRCA2" in df["approved_symbol"].values

    def test_fetch_acmg_incidental_data_disabled(self):
        """Test when ACMG is disabled."""
        config = {"data_sources": {"ACMG_Incidental_Findings": {"enabled": False}}}

        df = fetch_acmg_incidental_data(config)
        assert df.empty

    def test_fetch_acmg_genes_from_text_file(self):
        """Test loading ACMG genes from text file via fetch function."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\n")
            f.write("BRCA2\n")
            f.write("TP53\n")
            f.write("# This is a comment\n")
            f.write("\n")  # Empty line
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "ACMG_Incidental_Findings": {
                            "enabled": True,
                            "file_path": f.name,
                            "evidence_score": 1.0,
                        }
                    }
                }

                df = fetch_acmg_incidental_data(config)

                assert len(df) == 3
                assert "BRCA1" in df["approved_symbol"].values
                assert "BRCA2" in df["approved_symbol"].values
                assert "TP53" in df["approved_symbol"].values
            finally:
                Path(f.name).unlink()

    def test_fetch_acmg_genes_from_csv_file(self):
        """Test loading ACMG genes from CSV file via fetch function."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("Gene Symbol,Description\n")
            f.write("BRCA1,Breast cancer gene\n")
            f.write("BRCA2,Breast cancer gene\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "ACMG_Incidental_Findings": {
                            "enabled": True,
                            "file_path": f.name,
                            "evidence_score": 1.0,
                        }
                    }
                }

                df = fetch_acmg_incidental_data(config)

                assert len(df) == 2
                assert "BRCA1" in df["approved_symbol"].values
                assert "BRCA2" in df["approved_symbol"].values
            finally:
                Path(f.name).unlink()

    def test_fetch_acmg_genes_nonexistent_file(self):
        """Test loading from nonexistent file falls back to default."""
        config = {
            "data_sources": {
                "acmg_incidental": {
                    "enabled": True,
                    "file_path": "/nonexistent/file.txt",
                    "evidence_score": 1.0,
                }
            }
        }

        df = fetch_acmg_incidental_data(config)
        # Should fall back to default list
        assert len(df) == len(DEFAULT_ACMG_GENES)

    def test_fetch_acmg_with_custom_file(self):
        """Test fetching ACMG data with custom file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CUSTOM1\nCUSTOM2\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "ACMG_Incidental_Findings": {
                            "enabled": True,
                            "file_path": f.name,
                            "evidence_score": 2.0,
                        }
                    }
                }

                df = fetch_acmg_incidental_data(config)

                assert len(df) == 2
                assert "CUSTOM1" in df["approved_symbol"].values
                assert "CUSTOM2" in df["approved_symbol"].values
                assert all(df["source_evidence_score"] == 2.0)
            finally:
                Path(f.name).unlink()

    def test_process_inhouse_panel_missing_gene_column(self):
        """Test process_inhouse_panel with missing gene_column."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create Excel file without expected gene column
            excel_path = Path(tmpdir) / "test.xlsx"
            df = pd.DataFrame({"wrong_column": ["BRCA1", "TP53"], "other_col": [1, 2]})
            df.to_excel(excel_path, index=False)

            config = {
                "file_path": str(excel_path),
                "format": "excel",
                "gene_column": "gene_symbol",  # This column doesn't exist
                "evidence_score": 1.0,
            }

            result = process_inhouse_panel(config)
            # Should return None when gene column is missing
            assert result is None

    def test_extract_genes_from_column_mixed_types(self):
        """Test extract_genes_from_column with mixed data types."""
        # Create DataFrame with mixed data types in gene column
        mixed_df = pd.DataFrame(
            {
                "gene_symbol": [
                    "BRCA1",  # Valid gene
                    123,  # Number
                    "TP53",  # Valid gene
                    None,  # NaN
                    "",  # Empty string
                    "EGFR_VARIANT",  # Gene with underscore
                    3.14,  # Float
                ]
            }
        )

        genes = extract_genes_from_column(mixed_df, "gene_symbol")

        # Should extract all non-null values converted to strings
        # Note: Current implementation converts all non-null values to strings
        expected_genes = ["BRCA1", "TP53", "EGFR_VARIANT", "123", "3.14"]
        assert sorted(genes) == sorted(expected_genes)

    def test_fetch_acmg_incidental_data_empty_file(self):
        """Test ACMG data fetching with empty file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            # Write empty file
            f.write("")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "acmg_incidental": {
                            "enabled": True,
                            "file_path": f.name,
                            "evidence_score": 1.5,
                        }
                    }
                }

                df = fetch_acmg_incidental_data(config)

                # Should fall back to default ACMG genes
                assert len(df) == len(DEFAULT_ACMG_GENES)
                assert all(gene in DEFAULT_ACMG_GENES for gene in df["approved_symbol"])
            finally:
                Path(f.name).unlink()

    def test_fetch_acmg_incidental_data_duplicate_genes(self):
        """Test ACMG data fetching with duplicate genes in custom file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            # Write file with duplicates
            f.write("BRCA1\nTP53\nBRCA1\nTP53\nEGFR\nBRCA1\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "ACMG_Incidental_Findings": {
                            "enabled": True,
                            "file_path": f.name,
                            "evidence_score": 1.5,
                        }
                    }
                }

                df = fetch_acmg_incidental_data(config)

                # Should contain only unique genes
                unique_genes = df["approved_symbol"].unique()
                assert len(unique_genes) == 3  # BRCA1, TP53, EGFR
                assert sorted(unique_genes) == ["BRCA1", "EGFR", "TP53"]
            finally:
                Path(f.name).unlink()

    def test_fetch_acmg_data_from_live_source(self):
        """Test fetching ACMG data from live URL with mocked response."""
        # Load mock HTML fixture
        fixture_path = (
            Path(__file__).parent / "data" / "scraper_fixtures" / "acmg_ncbi.html"
        )
        with open(fixture_path) as f:
            mock_html = f.read()

        # Mock the requests.get call
        mock_response = Mock()
        mock_response.content = mock_html.encode("utf-8")
        mock_response.raise_for_status.return_value = None

        config = {
            "data_sources": {
                "ACMG_Incidental_Findings": {
                    "enabled": True,
                    "url": "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/",
                    "evidence_score": 1.5,
                }
            }
        }

        with patch(
            "custom_panel.sources.a_incidental_findings.requests.get",
            return_value=mock_response,
        ):
            df = fetch_acmg_incidental_data(config)

        # The function may succeed with mock (5 genes), use default list (71 genes),
        # or scrape real data (variable count). Test should be flexible.

        # Verify we got some genes
        assert len(df) > 0, "Should return at least some genes"

        # Check that essential cancer genes are present (these should be in any ACMG list)
        essential_genes = ["BRCA1", "BRCA2", "TP53", "APC"]
        for gene in essential_genes:
            assert gene in df["approved_symbol"].values, (
                f"Essential gene {gene} missing"
            )

        # Verify data structure
        expected_columns = [
            "approved_symbol",
            "source_evidence_score",
            "source_name",
            "source_details",
        ]
        for col in expected_columns:
            assert col in df.columns, f"Missing expected column: {col}"

        # Check source details indicate URL source (whether mock or real)
        assert all("URL:" in detail for detail in df["source_details"].values)
        assert all(
            "ncbi.nlm.nih.gov" in detail for detail in df["source_details"].values
        )

    def test_fetch_acmg_data_scraper_fails(self):
        """Test ACMG data fetching when scraper fails, should fall back to default list."""
        config = {
            "data_sources": {
                "acmg_incidental": {
                    "enabled": True,
                    "url": "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/",
                    "evidence_score": 1.5,
                }
            }
        }

        # Mock requests.get to raise an exception
        with patch(
            "custom_panel.sources.a_incidental_findings.requests.get",
            side_effect=Exception("Network error"),
        ):
            df = fetch_acmg_incidental_data(config)

        # Should fall back to default list
        assert len(df) == len(DEFAULT_ACMG_GENES)
        assert all(gene in DEFAULT_ACMG_GENES for gene in df["approved_symbol"])

        # Check source details indicate default source
        assert all(
            "ACMG_SF_v3.2_default" in detail for detail in df["source_details"].values
        )

    def test_scrape_acmg_genes_from_ncbi_success(self):
        """Test successful scraping of ACMG genes from NCBI."""
        # Load mock HTML fixture
        fixture_path = (
            Path(__file__).parent / "data" / "scraper_fixtures" / "acmg_ncbi.html"
        )
        with open(fixture_path) as f:
            mock_html = f.read()

        # Mock the requests.get call
        mock_response = Mock()
        mock_response.content = mock_html.encode("utf-8")
        mock_response.raise_for_status.return_value = None

        with patch(
            "custom_panel.sources.a_incidental_findings.requests.get",
            return_value=mock_response,
        ):
            genes = _scrape_acmg_genes_from_ncbi(
                "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"
            )

        expected_genes = ["APC", "BRCA1", "BRCA2", "TP53", "PTEN"]
        assert genes == expected_genes

    def test_scrape_acmg_genes_from_ncbi_network_error(self):
        """Test scraping when network request fails."""
        with patch(
            "custom_panel.sources.a_incidental_findings.requests.get",
            side_effect=Exception("Network error"),
        ):
            try:
                _scrape_acmg_genes_from_ncbi(
                    "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"
                )
                raise AssertionError("Should have raised an exception")
            except Exception as e:
                assert "Network error" in str(e)

    def test_scrape_acmg_genes_from_ncbi_missing_table(self):
        """Test scraping when expected table structure is missing."""
        mock_html = (
            "<html><body><div id='maincontent'><p>No table here</p></div></body></html>"
        )

        mock_response = Mock()
        mock_response.content = mock_html.encode("utf-8")
        mock_response.raise_for_status.return_value = None

        with patch(
            "custom_panel.sources.a_incidental_findings.requests.get",
            return_value=mock_response,
        ):
            try:
                _scrape_acmg_genes_from_ncbi(
                    "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"
                )
                raise AssertionError("Should have raised ValueError")
            except ValueError as e:
                assert "Could not find the ACMG gene table" in str(e)


class TestManualCuration:
    """Test manual curation functionality."""

    def test_validate_manual_curation_config_valid(self):
        """Test validation of valid manual curation config."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\nTP53\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "manual_curation": {
                            "enabled": True,
                            "lists": [
                                {
                                    "name": "Test List",
                                    "file_path": f.name,
                                    "gene_column": "gene_symbol",
                                    "evidence_score": 1.0,
                                }
                            ],
                        }
                    }
                }

                errors = validate_manual_curation_config(config)
                assert len(errors) == 0
            finally:
                Path(f.name).unlink()

    def test_validate_manual_curation_config_invalid(self):
        """Test validation of invalid manual curation config."""
        config = {
            "data_sources": {
                "Manual_Curation": {
                    "enabled": True,
                    "lists": [
                        {
                            "file_path": "/nonexistent/file.txt",
                            "evidence_score": -1.0,
                        }
                    ],
                }
            }
        }

        errors = validate_manual_curation_config(config)
        assert len(errors) > 0
        assert any("name is required" in error for error in errors)
        assert any("does not exist" in error for error in errors)
        assert any("must be a number between 0 and 1" in error for error in errors)

    def test_process_manual_list_csv(self):
        """Test processing CSV manual list."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("Gene Symbol,Description\n")
            f.write("BRCA1,Breast cancer gene\n")
            f.write("TP53,Tumor suppressor\n")
            f.flush()

            try:
                list_config = {
                    "name": "Test_CSV_List",
                    "file_path": f.name,
                    "gene_column": "Gene Symbol",
                    "evidence_score": 0.9,
                }

                df = process_manual_list(list_config)

                assert df is not None
                assert len(df) == 2
                assert "BRCA1" in df["approved_symbol"].values
                assert "TP53" in df["approved_symbol"].values
                assert all(df["source_evidence_score"] == 0.9)
                assert all(df["source_name"] == "Manual_Curation:Test_CSV_List")
                assert "File:" in df.iloc[0]["source_details"]
                assert "Column:Gene Symbol" in df.iloc[0]["source_details"]
                # Test default category behavior
                assert "category" in df.columns
                assert all(df["category"] == "germline")  # Should default to germline
                assert "Category:germline" in df.iloc[0]["source_details"]
            finally:
                Path(f.name).unlink()

    def test_process_manual_list_excel(self):
        """Test processing Excel manual list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            excel_path = Path(tmpdir) / "test.xlsx"
            test_df = pd.DataFrame(
                {
                    "Gene": ["BRCA1", "TP53", "EGFR"],
                    "Description": ["Gene 1", "Gene 2", "Gene 3"],
                }
            )
            test_df.to_excel(excel_path, index=False)

            list_config = {
                "name": "Test_Excel_List",
                "file_path": str(excel_path),
                "gene_column": "Gene",
                "evidence_score": 1.0,
            }

            df = process_manual_list(list_config)

            assert df is not None
            assert len(df) == 3
            assert set(df["approved_symbol"]) == {"BRCA1", "TP53", "EGFR"}
            assert all(df["source_evidence_score"] == 1.0)
            assert all(df["source_name"] == "Manual_Curation:Test_Excel_List")

    def test_process_manual_list_txt(self):
        """Test processing text manual list."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\n")
            f.write("TP53\n")
            f.write("EGFR\n")
            f.flush()

            try:
                list_config = {
                    "name": "Test_Text_List",
                    "file_path": f.name,
                    "gene_column": "Gene",
                    "evidence_score": 0.8,
                }

                df = process_manual_list(list_config)

                assert df is not None
                assert len(df) == 3
                assert set(df["approved_symbol"]) == {"BRCA1", "TP53", "EGFR"}
                assert all(df["source_evidence_score"] == 0.8)
                assert all(df["source_name"] == "Manual_Curation:Test_Text_List")
            finally:
                Path(f.name).unlink()

    def test_process_manual_list_missing_file(self):
        """Test processing manual list with missing file."""
        list_config = {
            "name": "Missing_File_List",
            "file_path": "/nonexistent/file.txt",
            "gene_column": "Gene",
            "evidence_score": 1.0,
        }

        df = process_manual_list(list_config)
        assert df.empty

    def test_process_manual_list_unsupported_format(self):
        """Test processing manual list with unsupported format."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            f.write('{"genes": ["BRCA1", "TP53"]}')
            f.flush()

            try:
                list_config = {
                    "name": "JSON_List",
                    "file_path": f.name,
                    "gene_column": "Gene",
                    "evidence_score": 1.0,
                }

                df = process_manual_list(list_config)
                assert df.empty
            finally:
                Path(f.name).unlink()

    def test_fetch_manual_curation_data_multiple_lists(self):
        """Test fetching manual curation data with multiple lists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create first list file
            txt_path = Path(tmpdir) / "list1.txt"
            txt_path.write_text("BRCA1\nTP53\n")

            # Create second list file
            csv_path = Path(tmpdir) / "list2.csv"
            csv_df = pd.DataFrame({"Genes": ["EGFR", "KRAS"]})
            csv_df.to_csv(csv_path, index=False)

            config = {
                "data_sources": {
                    "Manual_Curation": {
                        "enabled": True,
                        "lists": [
                            {
                                "name": "List1",
                                "file_path": str(txt_path),
                                "gene_column": "Gene",
                                "evidence_score": 1.0,
                            },
                            {
                                "name": "List2",
                                "file_path": str(csv_path),
                                "gene_column": "Genes",
                                "evidence_score": 0.8,
                            },
                        ],
                    }
                }
            }

            df = fetch_manual_curation_data(config)

            assert len(df) == 4
            assert set(df["approved_symbol"]) == {"BRCA1", "TP53", "EGFR", "KRAS"}

            # Check evidence scores are preserved
            list1_genes = df[df["source_name"] == "Manual_Curation:List1"]
            list2_genes = df[df["source_name"] == "Manual_Curation:List2"]
            assert all(list1_genes["source_evidence_score"] == 1.0)
            assert all(list2_genes["source_evidence_score"] == 0.8)

    def test_fetch_manual_curation_data_disabled(self):
        """Test when manual curation is disabled."""
        config = {"data_sources": {"manual_curation": {"enabled": False}}}

        df = fetch_manual_curation_data(config)
        assert df.empty

    def test_fetch_manual_curation_data_no_lists(self):
        """Test when no lists are configured."""
        config = {"data_sources": {"manual_curation": {"enabled": True, "lists": []}}}

        df = fetch_manual_curation_data(config)
        assert df.empty

    def test_fetch_manual_curation_data_all_lists_fail(self):
        """Test when all configured lists fail to process."""
        config = {
            "data_sources": {
                "manual_curation": {
                    "enabled": True,
                    "lists": [
                        {
                            "name": "Missing1",
                            "file_path": "/nonexistent1.txt",
                            "gene_column": "Gene",
                            "evidence_score": 1.0,
                        },
                        {
                            "name": "Missing2",
                            "file_path": "/nonexistent2.txt",
                            "gene_column": "Gene",
                            "evidence_score": 1.0,
                        },
                    ],
                }
            }
        }

        df = fetch_manual_curation_data(config)
        assert df.empty

    def test_get_manual_curation_summary(self):
        """Test getting manual curation summary."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\nTP53\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "Manual_Curation": {
                            "enabled": True,
                            "lists": [
                                {
                                    "name": "Test List",
                                    "file_path": f.name,
                                    "gene_column": "gene_symbol",
                                    "evidence_score": 1.0,
                                }
                            ],
                        }
                    }
                }

                summary = get_manual_curation_summary(config)

                assert summary["enabled"] is True
                assert summary["lists_count"] == 1
                assert len(summary["validation_errors"]) == 0
                assert summary["total_genes"] == 2
            finally:
                Path(f.name).unlink()

    def test_process_manual_list_empty_file(self):
        """Test processing manual list with empty file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("")  # Empty file
            f.flush()

            try:
                list_config = {
                    "name": "Empty_List",
                    "file_path": f.name,
                    "gene_column": "Gene",
                    "evidence_score": 1.0,
                }

                df = process_manual_list(list_config)
                assert df.empty
            finally:
                Path(f.name).unlink()

    def test_process_manual_list_missing_gene_column(self):
        """Test processing manual list with missing gene column."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("Other_Column,Description\n")
            f.write("value1,desc1\n")
            f.write("value2,desc2\n")
            f.flush()

            try:
                list_config = {
                    "name": "Missing_Column_List",
                    "file_path": f.name,
                    "gene_column": "Gene_Symbol",  # This column doesn't exist
                    "evidence_score": 1.0,
                }

                df = process_manual_list(list_config)
                assert df.empty
            finally:
                Path(f.name).unlink()

    def test_process_manual_list_default_parameters(self):
        """Test processing manual list with default parameters."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("BRCA1\nTP53\n")
            f.flush()

            try:
                # Minimal config - should use defaults
                list_config = {
                    "name": "Default_List",
                    "file_path": f.name,
                }

                df = process_manual_list(list_config)

                assert df is not None
                assert len(df) == 2
                assert all(df["source_evidence_score"] == 1.0)  # Default evidence score
                assert (
                    "Column:gene_symbol" in df.iloc[0]["source_details"]
                )  # Default gene column
            finally:
                Path(f.name).unlink()


class TestCOSMICGermlineData:
    """Test COSMIC Cancer Gene Census germline functionality."""

    def test_calculate_cosmic_score(self):
        """Test COSMIC tier score calculation."""
        tier_weights = {"Tier 1": 1.0, "Tier 2": 0.8, "": 0.4}

        assert _calculate_cosmic_score("Tier 1", tier_weights) == 1.0
        assert _calculate_cosmic_score("Tier 2", tier_weights) == 0.8
        assert _calculate_cosmic_score("", tier_weights) == 0.4
        assert _calculate_cosmic_score("Unknown", tier_weights) == 0.4
        assert _calculate_cosmic_score(None, tier_weights) == 0.4

    def test_is_cache_valid(self):
        """Test cache validity checking."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            cache_path = Path(f.name)

        try:
            # Fresh file should be valid
            assert _is_cache_valid(cache_path, 30) is True

            # Non-existent file should be invalid
            cache_path.unlink()
            assert _is_cache_valid(cache_path, 30) is False

        finally:
            if cache_path.exists():
                cache_path.unlink()

    def test_validate_cosmic_config_valid(self):
        """Test validation of valid COSMIC germline config."""
        config = {
            "data_sources": {
                "cosmic_germline": {
                    "enabled": True,
                    "email": "test@example.com",
                    "password": "testpass",
                    "cache_expiry_days": 30,
                    "germline_scoring": {
                        "enabled": True,
                        "tier_weights": {"Tier 1": 1.0, "Tier 2": 0.8, "": 0.4},
                    },
                }
            }
        }

        errors = validate_cosmic_config(config)
        assert len(errors) == 0

    def test_validate_cosmic_config_invalid(self):
        """Test validation of invalid COSMIC germline config."""
        config = {
            "data_sources": {
                "COSMIC_Germline": {
                    "enabled": True,
                    "cache_expiry_days": -1,
                    "germline_scoring": {
                        "enabled": True,
                        "tier_weights": {"Tier 1": 2.0},  # Invalid weight > 1
                    },
                }
            }
        }

        errors = validate_cosmic_config(config)
        assert len(errors) > 0
        assert any("requires either credentials" in error for error in errors)
        assert any(
            "cache_expiry_days must be a positive integer" in error for error in errors
        )
        assert any(
            "tier weight" in error and "must be between 0 and 1" in error
            for error in errors
        )

    def test_fetch_cosmic_germline_data_disabled(self):
        """Test when COSMIC germline is disabled."""
        config = {"data_sources": {"cosmic_germline": {"enabled": False}}}

        df = fetch_cosmic_germline_data(config)
        assert df.empty

    @patch("custom_panel.sources.g04_cosmic_germline._download_cosmic_census")
    @patch("custom_panel.sources.g04_cosmic_germline._load_cosmic_census")
    def test_fetch_cosmic_germline_data_success(self, mock_load, mock_download):
        """Test successful COSMIC germline data fetching."""
        # Mock COSMIC census data
        mock_df = pd.DataFrame(
            {
                "Gene Symbol": ["BRCA1", "TP53", "ATM"],
                "Tier": ["Tier 1", "Tier 1", "Tier 2"],
                "Germline": ["yes", "yes", "yes"],
            }
        )
        mock_df.attrs = {"germline_col": "Germline"}
        mock_load.return_value = mock_df

        config = {
            "data_sources": {
                "COSMIC_Germline": {
                    "enabled": True,
                    "census_url": "https://example.com/census.csv",
                    "cache_dir": ".cache/cosmic",
                    "cache_expiry_days": 30,
                    "germline_scoring": {
                        "enabled": True,
                        "tier_weights": {"Tier 1": 1.0, "Tier 2": 0.8, "": 0.4},
                    },
                }
            }
        }

        df = fetch_cosmic_germline_data(config)

        # Should have germline entries only
        assert len(df) > 0
        assert all(df["category"] == "germline")
        assert "BRCA1" in df["approved_symbol"].values
        assert "TP53" in df["approved_symbol"].values
        assert "ATM" in df["approved_symbol"].values

    def test_get_cosmic_summary(self):
        """Test COSMIC germline configuration summary."""
        config = {
            "data_sources": {
                "COSMIC_Germline": {
                    "enabled": True,
                    "email": "test@example.com",
                    "password": "testpass",
                    "census_url": "https://example.com/census.csv",
                    "cache_dir": ".cache/cosmic_test",
                    "cache_expiry_days": 30,
                    "germline_scoring": {"enabled": True},
                }
            }
        }

        summary = get_cosmic_summary(config)

        assert summary["enabled"] is True
        assert summary["has_credentials"] is True
        assert summary["census_url"] == "https://example.com/census.csv"
        assert summary["cache_dir"] == ".cache/cosmic_test"
        assert summary["cache_expiry_days"] == 30
        assert summary["germline_enabled"] is True
        assert summary["cache_exists"] is False
        assert summary["cache_valid"] is False
