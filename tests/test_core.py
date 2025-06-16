"""
Tests for core modules (HGNC client, Ensembl client, I/O operations).
"""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from custom_panel.core.ensembl_client import EnsemblClient
from custom_panel.core.hgnc_client import HGNCClient
from custom_panel.core.io import (
    STANDARD_COLUMNS,
    create_bed_file,
    create_standard_dataframe,
    load_panel_data,
    save_panel_data,
    validate_panel_dataframe,
)


class TestHGNCClient:
    """Test HGNC client functionality."""

    def test_init(self):
        """Test client initialization."""
        client = HGNCClient(timeout=10, max_retries=2)
        assert client.timeout == 10
        assert client.max_retries == 2

    @patch("custom_panel.core.hgnc_client.requests.Session.post")
    def test_standardize_symbols_batch(self, mock_post):
        """Test batch symbol standardization."""
        # Mock successful batch response
        mock_response = Mock()
        mock_response.json.return_value = {
            "response": {
                "docs": [
                    {"symbol": "BRCA1", "alias_symbol": ["brca1"], "prev_symbol": []},
                    {"symbol": "TP53", "alias_symbol": ["tp53"], "prev_symbol": []},
                ]
            }
        }
        mock_response.raise_for_status.return_value = None
        mock_post.return_value = mock_response

        client = HGNCClient()
        result = client.standardize_symbols_batch(("brca1", "tp53", "notfound"))

        expected = {
            "brca1": "BRCA1",  # Maps to standardized symbol
            "tp53": "TP53",  # Maps to standardized symbol
            "notfound": "notfound",  # Not found, returns original (lowercase preserved)
        }
        assert result == expected
        mock_post.assert_called_once()

    @patch("custom_panel.core.hgnc_client.requests.Session.post")
    def test_standardize_symbols_batch_fallback(self, mock_post):
        """Test batch standardization with fallback to individual requests."""
        import requests

        # Mock failed batch response
        mock_post.side_effect = requests.RequestException("API error")

        with patch.object(HGNCClient, "standardize_symbol") as mock_individual:
            mock_individual.side_effect = lambda x: x.upper()

            client = HGNCClient()
            result = client.standardize_symbols_batch(("brca1", "tp53"))

            expected = {"brca1": "BRCA1", "tp53": "TP53"}
            assert result == expected
            assert mock_individual.call_count == 2

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_symbol_to_hgnc_id(self, mock_get):
        """Test symbol to HGNC ID conversion."""
        # Mock successful response
        mock_response = Mock()
        mock_response.json.return_value = {
            "response": {"docs": [{"hgnc_id": "HGNC:1100"}]}
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HGNCClient()
        result = client.symbol_to_hgnc_id("BRCA1")

        assert result == "HGNC:1100"
        mock_get.assert_called_once()

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_symbol_to_hgnc_id_request_exception(self, mock_get):
        """Test symbol lookup with request exception."""
        import requests

        mock_get.side_effect = requests.RequestException("Network error")

        client = HGNCClient()
        result = client.symbol_to_hgnc_id("BRCA1")

        assert result is None

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_symbol_to_hgnc_id_404_error(self, mock_get):
        """Test symbol lookup with 404 error."""
        mock_response = Mock()
        import requests

        mock_response.raise_for_status.side_effect = requests.RequestException(
            "404 Not Found"
        )
        mock_get.return_value = mock_response

        client = HGNCClient()
        result = client.symbol_to_hgnc_id("BRCA1")

        assert result is None

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_symbol_to_hgnc_id_malformed_json(self, mock_get):
        """Test symbol lookup with malformed JSON response."""
        mock_response = Mock()
        import json

        mock_response.json.side_effect = json.JSONDecodeError("Invalid JSON", "", 0)
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HGNCClient()
        result = client.symbol_to_hgnc_id("BRCA1")

        assert result is None

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_symbol_to_hgnc_id_missing_keys(self, mock_get):
        """Test symbol lookup with missing expected keys in response."""
        mock_response = Mock()
        mock_response.json.return_value = {"unexpected": "structure"}
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HGNCClient()
        result = client.symbol_to_hgnc_id("BRCA1")

        assert result is None

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_symbol_to_hgnc_id_not_found(self, mock_get):
        """Test symbol to HGNC ID when gene not found."""
        mock_response = Mock()
        mock_response.json.return_value = {"response": {"docs": []}}
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HGNCClient()
        result = client.symbol_to_hgnc_id("NONEXISTENT")

        assert result is None

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_standardize_symbol(self, mock_get):
        """Test gene symbol standardization."""
        # Mock response for direct lookup
        mock_response = Mock()
        mock_response.json.return_value = {"response": {"docs": [{"symbol": "BRCA1"}]}}
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = HGNCClient()
        result = client.standardize_symbol("brca1")

        assert result == "BRCA1"


class TestEnsemblClient:
    """Test Ensembl client functionality."""

    def test_init(self):
        """Test client initialization."""
        client = EnsemblClient(timeout=15, max_retries=3)
        assert client.timeout == 15
        assert client.max_retries == 3

    @patch("custom_panel.core.ensembl_client.requests.Session.get")
    def test_get_gene_coordinates(self, mock_get):
        """Test getting gene coordinates."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "id": "ENSG00000012048",
            "seq_region_name": "17",
            "start": 43044295,
            "end": 43125364,
            "strand": -1,
            "biotype": "protein_coding",
            "description": "BRCA1 DNA repair associated",
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = EnsemblClient()
        result = client.get_gene_coordinates("BRCA1")

        assert result is not None
        assert result["gene_id"] == "ENSG00000012048"
        assert result["chromosome"] == "17"
        assert result["start"] == 43044295
        assert result["end"] == 43125364

    @patch("custom_panel.core.ensembl_client.requests.Session.post")
    def test_get_genes_coordinates_batch(self, mock_post):
        """Test batch gene coordinates lookup."""
        mock_response = Mock()
        mock_response.json.return_value = {
            "BRCA1": {
                "id": "ENSG00000012048",
                "seq_region_name": "17",
                "start": 43044295,
                "end": 43125364,
                "strand": -1,
                "biotype": "protein_coding",
            },
            "TP53": {
                "id": "ENSG00000141510",
                "seq_region_name": "17",
                "start": 7661779,
                "end": 7687550,
                "strand": -1,
                "biotype": "protein_coding",
            },
        }
        mock_response.raise_for_status.return_value = None
        mock_post.return_value = mock_response

        client = EnsemblClient()
        result = client.get_genes_coordinates(["BRCA1", "TP53"])

        assert len(result) == 2
        assert "BRCA1" in result
        assert "TP53" in result
        assert result["BRCA1"]["gene_id"] == "ENSG00000012048"

    @patch("custom_panel.core.ensembl_client.requests.Session.get")
    def test_get_gene_coordinates_request_exception(self, mock_get):
        """Test gene coordinates lookup with request exception."""
        import requests

        mock_get.side_effect = requests.RequestException("Network error")

        client = EnsemblClient()
        result = client.get_gene_coordinates("BRCA1")

        assert result is None

    @patch("custom_panel.core.ensembl_client.requests.Session.get")
    def test_get_gene_coordinates_500_error(self, mock_get):
        """Test gene coordinates lookup with server error."""
        mock_response = Mock()
        import requests

        mock_response.raise_for_status.side_effect = requests.RequestException(
            "500 Server Error"
        )
        mock_get.return_value = mock_response

        client = EnsemblClient()
        result = client.get_gene_coordinates("BRCA1")

        assert result is None

    @patch("custom_panel.core.ensembl_client.requests.Session.get")
    def test_get_gene_coordinates_malformed_json(self, mock_get):
        """Test gene coordinates lookup with malformed JSON."""
        mock_response = Mock()
        import json

        mock_response.json.side_effect = json.JSONDecodeError("Invalid JSON", "", 0)
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = EnsemblClient()
        result = client.get_gene_coordinates("BRCA1")

        assert result is None

    @patch("custom_panel.core.ensembl_client.requests.Session.get")
    def test_get_gene_coordinates_missing_keys(self, mock_get):
        """Test gene coordinates lookup with missing expected keys."""
        mock_response = Mock()
        mock_response.json.return_value = {"unexpected": "structure"}
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        client = EnsemblClient()
        result = client.get_gene_coordinates("BRCA1")

        assert result is not None
        assert result["gene_id"] is None
        assert result["chromosome"] is None

    @patch("custom_panel.core.ensembl_client.requests.Session.get")
    @patch("custom_panel.core.ensembl_client.requests.Session.post")
    def test_get_genes_coordinates_batch_fails_fallback(self, mock_post, mock_get):
        """Test batch gene coordinates lookup with fallback to individual requests."""
        import requests

        # Mock batch POST request failure
        mock_post.side_effect = requests.RequestException("Batch API error")

        # Mock successful individual GET requests
        mock_get_response = Mock()
        mock_get_response.json.return_value = {
            "id": "ENSG00000012048",
            "seq_region_name": "17",
            "start": 43044295,
            "end": 43125364,
            "strand": -1,
            "biotype": "protein_coding",
        }
        mock_get_response.raise_for_status.return_value = None
        mock_get.return_value = mock_get_response

        client = EnsemblClient()
        result = client.get_genes_coordinates(["BRCA1", "TP53"])

        # Verify batch request was attempted (with retries: 1 initial + 3 retries = 4 total)
        assert mock_post.call_count == 4

        # Verify fallback to individual requests occurred
        assert mock_get.call_count == 2

        # Verify results structure
        assert len(result) == 2
        assert "BRCA1" in result
        assert "TP53" in result
        assert result["BRCA1"]["gene_id"] == "ENSG00000012048"
        assert result["TP53"]["gene_id"] == "ENSG00000012048"


class TestIO:
    """Test I/O operations."""

    def test_create_standard_dataframe(self):
        """Test creating standardized DataFrame."""
        genes = ["BRCA1", "TP53", "EGFR"]
        source_name = "Test_Panel"
        evidence_scores = [1.0, 0.8, 0.9]

        df = create_standard_dataframe(
            genes=genes, source_name=source_name, evidence_scores=evidence_scores
        )

        assert len(df) == 3
        assert list(df["approved_symbol"]) == genes
        assert all(df["source_name"] == source_name)
        assert list(df["source_evidence_score"]) == evidence_scores

        # Check all required columns are present
        for col in STANDARD_COLUMNS:
            assert col in df.columns

    def test_validate_panel_dataframe_valid(self):
        """Test DataFrame validation with valid data."""
        df = create_standard_dataframe(genes=["BRCA1", "TP53"], source_name="Test")

        # Should not raise exception
        assert validate_panel_dataframe(df) is True

    def test_validate_panel_dataframe_invalid(self):
        """Test DataFrame validation with invalid data."""
        # Empty DataFrame
        df = pd.DataFrame()
        with pytest.raises(ValueError, match="DataFrame is empty"):
            validate_panel_dataframe(df)

        # Missing required columns
        df = pd.DataFrame({"gene": ["BRCA1"]})
        with pytest.raises(ValueError, match="Missing required columns"):
            validate_panel_dataframe(df)

        # Invalid evidence scores
        df = create_standard_dataframe(["BRCA1"], "Test", [2.0])  # Score > 1.0
        with pytest.raises(ValueError, match="source_evidence_score must be between"):
            validate_panel_dataframe(df)

    def test_save_and_load_panel_data(self):
        """Test saving and loading panel data."""
        df = create_standard_dataframe(
            genes=["BRCA1", "TP53"], source_name="Test_Panel"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Test Parquet format
            parquet_path = Path(tmpdir) / "test.parquet"
            save_panel_data(df, parquet_path, "parquet")
            loaded_df = load_panel_data(parquet_path)

            pd.testing.assert_frame_equal(df, loaded_df)

            # Test CSV format
            csv_path = Path(tmpdir) / "test.csv"
            save_panel_data(df, csv_path, "csv")
            loaded_csv_df = load_panel_data(csv_path)

            # CSV might have different dtypes, so check content
            assert len(loaded_csv_df) == len(df)
            assert list(loaded_csv_df["approved_symbol"]) == list(df["approved_symbol"])

    def test_create_bed_file(self):
        """Test BED file creation."""
        # Create test data with coordinates
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "TP53", "EGFR"],
                "chromosome": ["17", "17", "7"],
                "gene_start": [43044295, 7661779, 55019021],
                "gene_end": [43125364, 7687550, 55211628],
                "gene_strand": [-1, -1, 1],
                "include_germline": [True, True, False],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "test.bed"
            create_bed_file(df, bed_path, "include_germline")

            # Read BED file and verify content
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 2  # Only 2 genes with include_germline=True

            # Check first line (TP53 - has smaller start position)
            fields = lines[0].strip().split("\t")
            assert fields[0] == "17"  # chromosome
            assert fields[1] == "7661778"  # start (0-based)
            assert fields[2] == "7687550"  # end
            assert fields[3] == "TP53"  # name

            # Check second line (BRCA1)
            fields = lines[1].strip().split("\t")
            assert fields[0] == "17"  # chromosome
            assert fields[1] == "43044294"  # start (0-based)
            assert fields[2] == "43125364"  # end
            assert fields[3] == "BRCA1"  # name

    def test_create_bed_file_missing_coordinates(self):
        """Test BED file creation with missing coordinates."""
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "UNKNOWN"],
                "chromosome": ["17", None],
                "gene_start": [43044295, None],
                "gene_end": [43125364, None],
                "include_germline": [True, True],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "test.bed"
            create_bed_file(df, bed_path, "include_germline")

            # Should only include genes with valid coordinates
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 1  # Only BRCA1 has coordinates

    def test_load_panel_data_nonexistent_file(self):
        """Test loading data from non-existent file."""
        with pytest.raises(FileNotFoundError):
            load_panel_data("/nonexistent/path/file.parquet")

    def test_load_panel_data_malformed_csv(self):
        """Test loading malformed CSV file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            malformed_csv = Path(tmpdir) / "malformed.csv"
            # Create CSV with incorrect delimiter and malformed structure
            with open(malformed_csv, "w") as f:
                f.write("col1;col2;col3\n")  # Wrong delimiter
                f.write("val1;val2;val3\n")
                f.write("incomplete;row\n")  # Missing column

            # Should handle gracefully (pandas is quite forgiving)
            df = load_panel_data(malformed_csv, "csv")
            assert not df.empty  # pandas reads what it can

    def test_load_panel_data_unsupported_format(self):
        """Test loading data with unsupported format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            unsupported_file = Path(tmpdir) / "test.txt"
            unsupported_file.write_text("some content")

            with pytest.raises(ValueError, match="Unsupported format"):
                load_panel_data(unsupported_file, "txt")

    def test_save_panel_data_invalid_dataframe(self):
        """Test saving invalid DataFrame."""
        # Create DataFrame missing required columns
        invalid_df = pd.DataFrame({"invalid_col": [1, 2, 3]})

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.parquet"

            with pytest.raises(ValueError, match="Missing required columns"):
                save_panel_data(invalid_df, output_path)

    def test_create_bed_file_natural_sorting(self):
        """Test BED file creation with natural chromosome sorting."""
        # Create DataFrame with chromosomes in non-alphabetical order
        df = pd.DataFrame(
            {
                "approved_symbol": ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
                "chromosome": ["chr10", "chrX", "chr2", "chr1", "chrY"],
                "gene_start": [1000, 2000, 3000, 4000, 5000],
                "gene_end": [2000, 3000, 4000, 5000, 6000],
                "include_all": [True, True, True, True, True],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "test_natural_sort.bed"
            create_bed_file(df, bed_path, "include_all")

            # Read the resulting BED file
            with open(bed_path) as f:
                lines = f.readlines()

            # Extract chromosomes from BED file (first column)
            bed_chromosomes = [line.split("\t")[0] for line in lines]

            # Should be in natural order: chr1, chr2, chr10, chrX, chrY
            expected_order = ["chr1", "chr2", "chr10", "chrX", "chrY"]
            assert bed_chromosomes == expected_order


class TestDataFrameOperations:
    """Test DataFrame manipulation operations."""

    def test_merge_panel_dataframes(self):
        """Test merging multiple panel DataFrames."""
        from custom_panel.core.io import merge_panel_dataframes

        df1 = create_standard_dataframe(["BRCA1", "TP53"], "Source1")
        df2 = create_standard_dataframe(["TP53", "EGFR"], "Source2")

        merged = merge_panel_dataframes([df1, df2])

        assert len(merged) == 4  # All records preserved
        assert len(merged["approved_symbol"].unique()) == 3  # 3 unique genes
        assert len(merged["source_name"].unique()) == 2  # 2 sources

    def test_get_unique_genes(self):
        """Test getting unique genes from DataFrame."""
        from custom_panel.core.io import get_unique_genes

        df = create_standard_dataframe(["BRCA1", "TP53", "BRCA1"], "Test")
        unique_genes = get_unique_genes(df)

        assert len(unique_genes) == 2
        assert "BRCA1" in unique_genes
        assert "TP53" in unique_genes
