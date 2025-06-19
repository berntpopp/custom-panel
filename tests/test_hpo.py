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
    _extract_genes_from_omim,
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

        # Current HPO API doesn't provide gene associations directly
        assert len(genes) == 0

    @patch("custom_panel.core.hpo_client.requests.Session.get")
    def test_get_descendant_terms(self, mock_get):
        """Test getting descendant terms."""

        # Mock responses for descendants endpoint (preferred) and fallback to individual terms
        def side_effect(url, **kwargs):
            mock_response = Mock()
            mock_response.raise_for_status.return_value = None

            if "descendants" in url and "HP:0002664" in url:
                # Mock descendants endpoint response
                mock_response.json.return_value = [
                    {"id": "HP:0012531"},
                    {"id": "HP:0012532"},
                    {"id": "HP:0012533"},
                ]
            elif "HP:0002664" in url:  # Root term
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
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("# OMIM genemap2.txt\n")
            f.write(
                "1\t123\t456\t1p36.33\t\t191170\tTP53\tTumor protein p53\tTP53\t7157\t\t\tLi-Fraumeni syndrome, 151623 (3)\t\n"
            )
            f.flush()

            try:
                config = {
                    "data_sources": {
                        "HPO_Neoplasm": {
                            "enabled": True,
                            "omim_genemap2_path": f.name,
                            "base_evidence_score": 0.7,
                            "cache_expiry_days": 30,
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
                "HPO_Neoplasm": {
                    "enabled": True,
                    "omim_genemap2_path": "/nonexistent/file.txt",
                    "base_evidence_score": 1.5,  # Invalid range
                    "cache_expiry_days": -1,  # Invalid range
                }
            }
        }

        errors = validate_hpo_neoplasm_config(config)
        assert len(errors) > 0
        assert any("does not exist" in error for error in errors)
        assert any("base_evidence_score must be between" in error for error in errors)
        assert any(
            "cache_expiry_days must be a positive integer" in error for error in errors
        )

    def test_extract_genes_from_omim(self):
        """Test extracting genes from OMIM data."""
        # Create mock OMIM DataFrame
        omim_df = pd.DataFrame(
            {
                "approved_symbol": ["TP53", "BRCA1", "KRAS"],
                "phenotypes": [
                    "Li-Fraumeni syndrome, 151623 (3), Autosomal dominant; Colorectal cancer, 114500 (3)",
                    "Breast-ovarian cancer, familial, 1, 604370 (3), Autosomal dominant",
                    "Cardiofaciocutaneous syndrome, 115150 (3); Noonan syndrome, 163950 (3)",
                ],
                "entrez_gene_id": ["7157", "672", "3845"],
                "mim_number": ["191170", "113705", "190070"],
            }
        )

        # Create mock phenotype DataFrame
        phenotype_df = pd.DataFrame(
            {
                "database_id": ["OMIM:151623", "OMIM:604370", "OMIM:115150"],
                "hpo_id": ["HP:0002664", "HP:0000729", "HP:0001631"],
                "disease_name": [
                    "Li-Fraumeni syndrome",
                    "Breast cancer",
                    "Cardiofaciocutaneous syndrome",
                ],
            }
        )

        config = {"base_evidence_score": 0.7}
        omim_ids = ["151623", "604370"]

        genes = _extract_genes_from_omim(omim_df, omim_ids, phenotype_df, config)

        assert len(genes) == 2
        assert "TP53" in genes
        assert "BRCA1" in genes
        assert "151623" in genes["TP53"]["omim_ids"]
        assert "604370" in genes["BRCA1"]["omim_ids"]

    def test_extract_genes_from_omim_no_matches(self):
        """Test extracting genes with no OMIM ID matches."""
        omim_df = pd.DataFrame(
            {
                "approved_symbol": ["TEST1"],
                "phenotypes": ["Test disease, 999999 (3)"],
            }
        )

        phenotype_df = pd.DataFrame(
            {
                "database_id": ["OMIM:888888"],
                "hpo_id": ["HP:0000001"],
                "disease_name": ["Different disease"],
            }
        )

        config = {"base_evidence_score": 0.7}
        omim_ids = ["151623"]

        genes = _extract_genes_from_omim(omim_df, omim_ids, phenotype_df, config)
        assert len(genes) == 0

    def test_extract_genes_from_omim_empty_phenotypes(self):
        """Test extracting genes with empty phenotypes."""
        omim_df = pd.DataFrame(
            {
                "approved_symbol": ["TEST1"],
                "phenotypes": [""],
                "entrez_gene_id": ["12345"],
                "mim_number": ["600000"],
            }
        )

        phenotype_df = pd.DataFrame(
            {
                "database_id": ["OMIM:151623"],
                "hpo_id": ["HP:0002664"],
                "disease_name": ["Some disease"],
            }
        )
        config = {"base_evidence_score": 0.7}
        omim_ids = ["151623"]

        genes = _extract_genes_from_omim(omim_df, omim_ids, phenotype_df, config)
        assert len(genes) == 0

    def test_calculate_evidence_score(self):
        """Test evidence score calculation."""
        config = {"base_evidence_score": 0.7}

        # Gene with multiple HPO terms and diseases
        gene_data = {
            "hpo_terms": {"HP:0002664", "HP:0012531", "HP:0001234"},  # 3 terms
            "diseases": {"Disease1", "Disease2", "Disease3"},  # 3 diseases
            "entrez_id": "7157",
        }

        score = _calculate_evidence_score(gene_data, config)

        # Base (0.7) + hpo bonus (3 * 0.02 = 0.06) + disease bonus (3 * 0.02 = 0.06) + entrez bonus (0.05) = 0.87
        expected_score = 0.7 + 0.06 + 0.06 + 0.05
        assert abs(score - expected_score) < 0.01

        # Gene with minimal data
        minimal_gene_data = {
            "hpo_terms": set(),
            "diseases": set(),
            "entrez_id": None,
        }

        minimal_score = _calculate_evidence_score(minimal_gene_data, config)
        # Just base score: 0.7
        assert minimal_score == 0.7

    def test_create_source_details(self):
        """Test creating source details string."""
        gene_data = {
            "hpo_terms": {"HP:0002664", "HP:0012531"},  # 2 terms
            "diseases": {"Disease1", "Disease2"},  # 2 diseases
            "omim_ids": {"151623", "604370"},  # 2 OMIM IDs
            "mim_number": "191170",
            "entrez_id": "7157",
        }

        details = _create_source_details(gene_data)

        assert "HPO_terms:2" in details
        assert "Diseases:2" in details
        assert "OMIM_IDs:2" in details
        assert "MIM:191170" in details
        assert "Entrez:7157" in details

    @patch("custom_panel.sources.g02_hpo.HPOClient")
    def test_fetch_hpo_neoplasm_data_disabled(self, mock_client_class):
        """Test when HPO neoplasm is disabled."""
        config = {"data_sources": {"hpo_neoplasm": {"enabled": False}}}

        df = fetch_hpo_neoplasm_data(config)
        assert df.empty

    @patch("custom_panel.sources.g02_hpo.HPOClient")
    def test_fetch_hpo_neoplasm_data_neoplasm_search(self, mock_client_class):
        """Test fetching data - should return empty since no OMIM genemap2 is configured."""
        # Mock HPO client
        mock_client = Mock()
        mock_client_class.return_value = mock_client

        # Mock descendant terms
        mock_client.get_descendant_terms.return_value = {"HP:0002664", "HP:0012531"}

        config = {
            "data_sources": {
                "hpo_neoplasm": {
                    "enabled": True,
                    "base_evidence_score": 0.5,
                }
            }
        }

        df = fetch_hpo_neoplasm_data(config)

        # Should return empty since no OMIM genemap2 URL or path is configured
        assert df.empty

    @patch("custom_panel.sources.g02_hpo.HPOClient")
    def test_fetch_hpo_neoplasm_data_specific_terms(self, mock_client_class):
        """Test fetching data - should return empty since no OMIM genemap2 is configured."""
        # Mock HPO client
        mock_client = Mock()
        mock_client_class.return_value = mock_client

        # Mock descendant terms
        mock_client.get_descendant_terms.return_value = {"HP:0002664"}

        config = {
            "data_sources": {
                "hpo_neoplasm": {
                    "enabled": True,
                    "base_evidence_score": 0.7,
                }
            }
        }

        df = fetch_hpo_neoplasm_data(config)

        # Should return empty since no OMIM genemap2 URL or path is configured
        assert df.empty

    def test_fetch_hpo_neoplasm_data_omim_file(self):
        """Test fetching data using OMIM genemap2 file."""
        # Create temporary OMIM genemap2 file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as omim_f:
            omim_f.write("# OMIM genemap2.txt\n")
            omim_f.write(
                "17\t43044295\t43125483\t17q21.31\t\t191170\tTP53\tTumor protein p53\tTP53\t7157\t\t\tLi-Fraumeni syndrome, 151623 (3); Colorectal cancer, 114500 (3)\t\n"
            )
            omim_f.write(
                "17\t43044295\t43125483\t17q21.31\t\t113705\tBRCA1\tBRCA1 DNA repair associated\tBRCA1\t672\t\t\tBreast-ovarian cancer, familial, 1, 604370 (3)\t\n"
            )
            omim_f.flush()

            # Create temporary phenotype.hpoa file
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".hpoa", delete=False
            ) as pheno_f:
                pheno_f.write(
                    "#database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration\n"
                )
                pheno_f.write(
                    "OMIM:151623\tLi-Fraumeni syndrome\t\tHP:0002664\tPMID:123\tTAS\t\t\t\t\tP\t\n"
                )
                pheno_f.write(
                    "OMIM:604370\tBreast-ovarian cancer, familial, 1\t\tHP:0002664\tPMID:456\tTAS\t\t\t\t\tP\t\n"
                )
                pheno_f.flush()

                try:
                    config = {
                        "data_sources": {
                            "HPO_Neoplasm": {
                                "enabled": True,
                                "omim_genemap2_path": omim_f.name,
                                "base_evidence_score": 0.6,
                                "cache_dir": Path(pheno_f.name).parent,
                            }
                        }
                    }

                    # Copy phenotype file to expected cache location
                    cache_dir = Path(pheno_f.name).parent
                    expected_pheno_path = cache_dir / "phenotype.hpoa"
                    import shutil

                    shutil.copy2(pheno_f.name, expected_pheno_path)

                    with patch(
                        "custom_panel.sources.g02_hpo.HPOClient"
                    ) as mock_client_class:
                        mock_client = Mock()
                        mock_client_class.return_value = mock_client

                        # Mock get_descendant_terms to return neoplasm terms
                        mock_client.get_descendant_terms.return_value = {"HP:0002664"}

                        # Mock file validity and parsing
                        mock_client.is_cache_valid.return_value = True

                        # Use real parsing methods
                        from custom_panel.core.hpo_client import HPOClient

                        mock_client.parse_phenotype_hpoa = (
                            HPOClient.parse_phenotype_hpoa
                        )
                        mock_client.parse_omim_genemap2 = HPOClient.parse_omim_genemap2

                        df = fetch_hpo_neoplasm_data(config)

                    assert len(df) == 2
                    assert "TP53" in df["approved_symbol"].values
                    assert "BRCA1" in df["approved_symbol"].values
                    assert all(df["source_name"] == "HPO_Neoplasm")

                finally:
                    Path(omim_f.name).unlink(missing_ok=True)
                    Path(pheno_f.name).unlink(missing_ok=True)
                    expected_pheno_path.unlink(missing_ok=True)

    def test_get_hpo_neoplasm_summary(self):
        """Test getting HPO neoplasm summary."""
        config = {
            "data_sources": {
                "HPO_Neoplasm": {
                    "enabled": True,
                    "omim_genemap2_path": "/path/to/omim.txt",
                    "base_evidence_score": 0.8,
                    "cache_expiry_days": 15,
                }
            }
        }

        summary = get_hpo_neoplasm_summary(config)

        assert summary["enabled"] is True
        assert summary["omim_genemap2_configured"] is True
        assert summary["base_evidence_score"] == 0.8
        assert summary["cache_expiry_days"] == 15
        assert "phenotype_hpoa_url" in summary
        assert "neoplasm_root_term" in summary

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
