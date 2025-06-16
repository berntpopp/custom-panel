"""
Tests for HPO client and HPO neoplasm functionality.
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd

from custom_panel.core.hpo_client import HPOClient
from custom_panel.sources.g02_hpo import (
    _calculate_evidence_score,
    _create_source_details,
    _get_genes_from_omim_file,
    fetch_hpo_neoplasm_data,
    get_hpo_neoplasm_summary,
    validate_hpo_neoplasm_config,
)


class TestHPOClient:
    """Test HPO client functionality."""

    def test_init(self):
        """Test client initialization."""
        client = HPOClient(timeout=10, max_retries=2, retry_delay=0.5)
        assert client.timeout == 10
        assert client.max_retries == 2
        assert client.retry_delay == 0.5

    @patch("custom_panel.core.hpo_client.requests.Session.get")
    def test_get_term_info(self, mock_get):
        """Test getting HPO term information."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "id": "HP:0002664",
            "name": "Neoplasm",
            "definition": "An organ or tissue overgrowth...",
            "synonyms": ["tumor", "cancer"],
            "parents": [{"id": "HP:0000001"}],
            "children": [{"id": "HP:0012531"}],
            "genes": [{"gene_symbol": "TP53"}],
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HPOClient()
        term_info = client.get_term_info("HP:0002664")

        assert term_info is not None
        assert term_info["id"] == "HP:0002664"
        assert term_info["name"] == "Neoplasm"
        assert "tumor" in term_info["synonyms"]

    @patch("custom_panel.core.hpo_client.requests.Session.get")
    def test_search_terms(self, mock_get):
        """Test searching for HPO terms."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "terms": [
                {
                    "id": "HP:0002664",
                    "name": "Neoplasm",
                    "definition": "An organ or tissue overgrowth...",
                    "synonyms": ["tumor", "cancer"],
                },
                {
                    "id": "HP:0012531",
                    "name": "Pain",
                    "definition": "An unpleasant sensation...",
                    "synonyms": [],
                },
            ]
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HPOClient()
        terms = client.search_terms("neoplasm")

        assert len(terms) == 2
        assert terms[0]["id"] == "HP:0002664"
        assert terms[0]["name"] == "Neoplasm"

    @patch("custom_panel.core.hpo_client.requests.Session.get")
    def test_get_term_genes(self, mock_get):
        """Test getting genes for an HPO term."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "genes": [
                {
                    "gene_symbol": "TP53",
                    "entrez_id": "7157",
                    "disease_id": "OMIM:151623",
                    "disease_name": "Li-Fraumeni syndrome",
                    "frequency": "HP:0040281",
                    "onset": "HP:0003577",
                    "evidence": "IEA",
                },
                {
                    "gene_symbol": "BRCA1",
                    "entrez_id": "672",
                    "disease_id": "OMIM:113705",
                    "disease_name": "Breast cancer",
                    "frequency": "HP:0040280",
                    "onset": "HP:0003581",
                    "evidence": "TAS",
                },
            ]
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HPOClient()
        genes = client.get_term_genes("HP:0002664")

        assert len(genes) == 2
        assert genes[0]["gene_symbol"] == "TP53"
        assert genes[0]["entrez_id"] == "7157"
        assert genes[1]["gene_symbol"] == "BRCA1"

    @patch("custom_panel.core.hpo_client.requests.Session.get")
    def test_get_descendant_terms(self, mock_get):
        """Test getting descendant terms."""

        # Mock responses for term hierarchy
        def side_effect(url, **kwargs):
            mock_response = Mock()
            mock_response.raise_for_status.return_value = None

            if "HP:0002664" in url:  # Root term
                mock_response.json.return_value = {
                    "id": "HP:0002664",
                    "name": "Neoplasm",
                    "children": [{"id": "HP:0012531"}, {"id": "HP:0012532"}],
                }
            elif "HP:0012531" in url:  # Child term 1
                mock_response.json.return_value = {
                    "id": "HP:0012531",
                    "name": "Child1",
                    "children": [{"id": "HP:0012533"}],
                }
            elif "HP:0012532" in url:  # Child term 2
                mock_response.json.return_value = {
                    "id": "HP:0012532",
                    "name": "Child2",
                    "children": [],
                }
            elif "HP:0012533" in url:  # Grandchild term
                mock_response.json.return_value = {
                    "id": "HP:0012533",
                    "name": "Grandchild",
                    "children": [],
                }

            return mock_response

        mock_get.side_effect = side_effect

        client = HPOClient()
        descendants = client.get_descendant_terms("HP:0002664", max_depth=3)

        # Should include the original term and all descendants
        expected_terms = {"HP:0002664", "HP:0012531", "HP:0012532", "HP:0012533"}
        assert descendants == expected_terms

    @patch("custom_panel.core.hpo_client.requests.Session.get")
    def test_find_neoplasm_terms(self, mock_get):
        """Test finding neoplasm-related terms."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "terms": [
                {
                    "id": "HP:0002664",
                    "name": "Neoplasm",
                    "definition": "An organ or tissue overgrowth...",
                    "synonyms": ["tumor"],
                },
                {
                    "id": "HP:0012531",
                    "name": "Tumor",
                    "definition": "A tumor is...",
                    "synonyms": ["neoplasm"],
                },
            ]
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HPOClient()
        terms = client.find_neoplasm_terms()

        # Should return unique terms (no duplicates)
        assert len(terms) >= 1
        term_ids = {term["id"] for term in terms}
        assert "HP:0002664" in term_ids or "HP:0012531" in term_ids

    def test_get_cache_info(self):
        """Test getting cache information."""
        client = HPOClient()
        cache_info = client.get_cache_info()

        assert "get_term_info" in cache_info
        assert "search_terms" in cache_info
        assert "get_term_genes" in cache_info
        assert "find_neoplasm_terms" in cache_info

    def test_clear_cache(self):
        """Test clearing cache."""
        client = HPOClient()
        # This should not raise an exception
        client.clear_cache()


class TestHPONeoplasm:
    """Test HPO neoplasm functionality."""

    def test_validate_hpo_neoplasm_config_valid(self):
        """Test validation of valid HPO neoplasm config."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("Gene Symbol,Disease\n")
            f.write("TP53,Li-Fraumeni syndrome\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "hpo_neoplasm": {
                            "enabled": True,
                            "omim_file_path": f.name,
                            "specific_hpo_terms": ["HP:0002664", "HP:0012531"],
                            "base_evidence_score": 0.7,
                            "max_hierarchy_depth": 5,
                        }
                    }
                }

                errors = validate_hpo_neoplasm_config(config)
                assert len(errors) == 0
            finally:
                Path(f.name).unlink()

    def test_validate_hpo_neoplasm_config_invalid(self):
        """Test validation of invalid HPO neoplasm config."""
        config = {
            "data_sources": {
                "hpo_neoplasm": {
                    "enabled": True,
                    "omim_file_path": "/nonexistent/file.csv",
                    "specific_hpo_terms": ["INVALID_ID", "HP:0002664"],
                    "base_evidence_score": 1.5,  # Invalid range
                    "max_hierarchy_depth": 25,  # Invalid range
                }
            }
        }

        errors = validate_hpo_neoplasm_config(config)
        assert len(errors) > 0
        assert any("does not exist" in error for error in errors)
        assert any("must be a valid HPO ID" in error for error in errors)
        assert any("base_evidence_score must be between" in error for error in errors)
        assert any("max_hierarchy_depth must be between" in error for error in errors)

    def test_get_genes_from_omim_file_csv(self):
        """Test getting genes from OMIM CSV file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("Gene Symbol,Disease,OMIM ID\n")
            f.write("TP53,Li-Fraumeni syndrome,151623\n")
            f.write("BRCA1,Breast cancer,113705\n")
            f.flush()

            try:
                config = {
                    "omim_gene_column": "Gene Symbol",
                    "omim_disease_column": "Disease",
                    "omim_id_column": "OMIM ID",
                }

                genes = _get_genes_from_omim_file(f.name, config)

                assert len(genes) == 2
                assert "TP53" in genes
                assert "BRCA1" in genes
                assert "Li-Fraumeni syndrome" in genes["TP53"]["diseases"]
                assert "OMIM" in genes["TP53"]["evidence_sources"]
                assert "omim_file" in genes["TP53"]["source_methods"]
            finally:
                Path(f.name).unlink()

    def test_get_genes_from_omim_file_excel(self):
        """Test getting genes from OMIM Excel file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            excel_path = Path(tmpdir) / "omim_data.xlsx"
            df = pd.DataFrame(
                {
                    "Gene": ["TP53", "BRCA1", "EGFR"],
                    "Disease": ["Li-Fraumeni", "Breast cancer", "Lung cancer"],
                    "OMIM": ["151623", "113705", "131550"],
                }
            )
            df.to_excel(excel_path, index=False)

            config = {
                "omim_gene_column": "Gene",
                "omim_disease_column": "Disease",
                "omim_id_column": "OMIM",
            }

            genes = _get_genes_from_omim_file(str(excel_path), config)

            assert len(genes) == 3
            assert "TP53" in genes
            assert "BRCA1" in genes
            assert "EGFR" in genes

    def test_get_genes_from_omim_file_missing(self):
        """Test getting genes from missing OMIM file."""
        config = {}
        genes = _get_genes_from_omim_file("/nonexistent/file.csv", config)
        assert len(genes) == 0

    def test_calculate_evidence_score(self):
        """Test evidence score calculation."""
        config = {"base_evidence_score": 0.5}

        # Gene with multiple data sources and HPO terms
        gene_data = {
            "source_methods": {"neoplasm_search", "specific_terms", "omim_file"},
            "hpo_terms": [{"hpo_id": "HP:0002664"}, {"hpo_id": "HP:0012531"}],
            "diseases": {"Disease1", "Disease2", "Disease3"},
            "entrez_id": "7157",
        }

        score = _calculate_evidence_score(gene_data, config)

        # Base (0.5) + source bonus (0.3) + hpo bonus (0.1) + disease bonus (0.06) + entrez bonus (0.05) = 1.01, capped at 1.0
        assert score == 1.0

        # Gene with minimal data
        minimal_gene_data = {
            "source_methods": {"neoplasm_search"},
            "hpo_terms": [],
            "diseases": set(),
            "entrez_id": None,
        }

        minimal_score = _calculate_evidence_score(minimal_gene_data, config)
        # Base (0.5) + source bonus (0.1) = 0.6
        assert minimal_score == 0.6

    def test_create_source_details(self):
        """Test creating source details string."""
        gene_data = {
            "source_methods": {"neoplasm_search", "omim_file"},
            "hpo_terms": [{"hpo_id": "HP:0002664"}, {"hpo_id": "HP:0012531"}],
            "diseases": {"Disease1", "Disease2"},
            "evidence_sources": {"HPO", "OMIM"},
            "entrez_id": "7157",
        }

        details = _create_source_details(gene_data)

        assert "Methods:neoplasm_search,omim_file" in details
        assert "HPO_terms:2" in details
        assert "Diseases:2" in details
        assert "Evidence:HPO,OMIM" in details
        assert "Entrez:7157" in details

    @patch("custom_panel.sources.g02_hpo.HPOClient")
    def test_fetch_hpo_neoplasm_data_disabled(self, mock_client_class):
        """Test when HPO neoplasm is disabled."""
        config = {"data_sources": {"hpo_neoplasm": {"enabled": False}}}

        df = fetch_hpo_neoplasm_data(config)
        assert df.empty

    @patch("custom_panel.sources.g02_hpo.HPOClient")
    def test_fetch_hpo_neoplasm_data_neoplasm_search(self, mock_client_class):
        """Test fetching data using neoplasm search."""
        # Mock HPO client
        mock_client = Mock()
        mock_client_class.return_value = mock_client

        # Mock neoplasm terms
        mock_client.find_neoplasm_terms.return_value = [
            {"id": "HP:0002664", "name": "Neoplasm"}
        ]

        # Mock hierarchy genes
        mock_client.get_genes_for_phenotype_hierarchy.return_value = {
            "TP53": {
                "gene_symbol": "TP53",
                "entrez_id": "7157",
                "hpo_terms": [{"hpo_id": "HP:0002664", "disease_name": "Li-Fraumeni"}],
                "diseases": ["Li-Fraumeni syndrome"],
                "evidence_sources": ["HPO"],
            }
        }

        config = {
            "data_sources": {
                "hpo_neoplasm": {
                    "enabled": True,
                    "use_neoplasm_search": True,
                    "base_evidence_score": 0.5,
                }
            }
        }

        df = fetch_hpo_neoplasm_data(config)

        assert len(df) == 1
        assert df.iloc[0]["approved_symbol"] == "TP53"
        assert df.iloc[0]["source_name"] == "HPO_OMIM_Neoplasm"
        assert df.iloc[0]["source_evidence_score"] >= 0.5

    @patch("custom_panel.sources.g02_hpo.HPOClient")
    def test_fetch_hpo_neoplasm_data_specific_terms(self, mock_client_class):
        """Test fetching data using specific HPO terms."""
        # Mock HPO client
        mock_client = Mock()
        mock_client_class.return_value = mock_client

        # Mock hierarchy genes for specific terms
        mock_client.get_genes_for_phenotype_hierarchy.return_value = {
            "BRCA1": {
                "gene_symbol": "BRCA1",
                "entrez_id": "672",
                "hpo_terms": [
                    {"hpo_id": "HP:0002664", "disease_name": "Breast cancer"}
                ],
                "diseases": ["Breast cancer"],
                "evidence_sources": ["HPO"],
            }
        }

        config = {
            "data_sources": {
                "hpo_neoplasm": {
                    "enabled": True,
                    "use_neoplasm_search": False,
                    "specific_hpo_terms": ["HP:0002664"],
                    "base_evidence_score": 0.7,
                }
            }
        }

        df = fetch_hpo_neoplasm_data(config)

        assert len(df) == 1
        assert df.iloc[0]["approved_symbol"] == "BRCA1"
        assert df.iloc[0]["source_evidence_score"] >= 0.7

    def test_fetch_hpo_neoplasm_data_omim_file(self):
        """Test fetching data using OMIM file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("Gene Symbol,Disease\n")
            f.write("TP53,Li-Fraumeni syndrome\n")
            f.write("BRCA1,Breast cancer\n")
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "hpo_neoplasm": {
                            "enabled": True,
                            "use_neoplasm_search": False,
                            "omim_file_path": f.name,
                            "base_evidence_score": 0.6,
                        }
                    }
                }

                with patch("custom_panel.sources.g02_hpo.HPOClient"):
                    df = fetch_hpo_neoplasm_data(config)

                assert len(df) == 2
                assert "TP53" in df["approved_symbol"].values
                assert "BRCA1" in df["approved_symbol"].values
                assert all(df["source_name"] == "HPO_OMIM_Neoplasm")
            finally:
                Path(f.name).unlink()

    def test_get_hpo_neoplasm_summary(self):
        """Test getting HPO neoplasm summary."""
        config = {
            "data_sources": {
                "hpo_neoplasm": {
                    "enabled": True,
                    "use_neoplasm_search": True,
                    "specific_hpo_terms": ["HP:0002664", "HP:0012531"],
                    "omim_file_path": "/path/to/omim.csv",
                }
            }
        }

        summary = get_hpo_neoplasm_summary(config)

        assert summary["enabled"] is True
        assert summary["use_neoplasm_search"] is True
        assert summary["specific_hpo_terms_count"] == 2
        assert summary["omim_file_configured"] is True
        assert "estimated_genes" in summary

    def test_fetch_hpo_neoplasm_data_no_sources(self):
        """Test when no data sources are configured."""
        config = {
            "data_sources": {
                "hpo_neoplasm": {
                    "enabled": True,
                    "use_neoplasm_search": False,
                    "specific_hpo_terms": [],
                    "omim_file_path": None,
                }
            }
        }

        with patch("custom_panel.sources.g02_hpo.HPOClient"):
            df = fetch_hpo_neoplasm_data(config)

        assert df.empty
