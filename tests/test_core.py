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
    create_exon_bed_file,
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

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_standardize_symbols_batch(self, mock_get):
        """Test batch symbol standardization."""
        # Mock successful batch response
        mock_response = Mock()
        mock_response.json.return_value = {
            "response": {
                "docs": [
                    {
                        "symbol": "BRCA1",
                        "hgnc_id": "HGNC:1100",
                        "alias_symbol": ["brca1"],
                        "prev_symbol": [],
                    },
                    {
                        "symbol": "TP53",
                        "hgnc_id": "HGNC:11998",
                        "alias_symbol": ["tp53"],
                        "prev_symbol": [],
                    },
                ]
            }
        }
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response

        # Need to also mock individual symbol lookups for "notfound"
        with patch.object(HGNCClient, "standardize_symbol") as mock_standardize:
            mock_standardize.return_value = "notfound"

            client = HGNCClient()
            result = client.standardize_symbols_batch(("brca1", "tp53", "notfound"))

            expected = {
                "brca1": {
                    "approved_symbol": "BRCA1",
                    "hgnc_id": "HGNC:1100",
                },  # Maps to standardized symbol with ID
                "tp53": {
                    "approved_symbol": "TP53",
                    "hgnc_id": "HGNC:11998",
                },  # Maps to standardized symbol with ID
                "notfound": {
                    "approved_symbol": "notfound",
                    "hgnc_id": None,
                },  # Not found, returns original (lowercase preserved)
            }
            assert result == expected
            mock_get.assert_called_once()
            mock_standardize.assert_called_once_with("notfound")

    @patch("custom_panel.core.hgnc_client.requests.Session.get")
    def test_standardize_symbols_batch_fallback(self, mock_get):
        """Test batch standardization with fallback to individual requests."""
        import requests

        # Mock failed batch response
        mock_get.side_effect = requests.RequestException("API error")

        with patch.object(HGNCClient, "standardize_symbol") as mock_individual:
            mock_individual.side_effect = lambda x: x.upper()

            with patch.object(HGNCClient, "get_gene_info") as mock_gene_info:
                # Mock gene info to return HGNC IDs
                mock_gene_info.side_effect = lambda x: {"hgnc_id": f"HGNC:{x}"}

                client = HGNCClient()
                result = client.standardize_symbols_batch(("brca1", "tp53"))

                expected = {
                    "brca1": {"approved_symbol": "BRCA1", "hgnc_id": "HGNC:BRCA1"},
                    "tp53": {"approved_symbol": "TP53", "hgnc_id": "HGNC:TP53"},
                }
                assert result == expected
                assert mock_individual.call_count == 2
                assert mock_gene_info.call_count == 2

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
        df = create_standard_dataframe(["BRCA1"], "Test", [6.0])  # Score > 5.0
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

    def test_create_exon_bed_file(self):
        """Test the creation of a BED file from exon data."""
        exons_data = [
            # Gene 1, canonical transcript
            {
                "chromosome": "1",
                "start": 100,
                "end": 200,
                "gene_symbol": "GENE1",
                "exon_id": "E1",
                "strand": "+",
                "transcript_id": "T1",
                "rank": 1,
            },
            {
                "chromosome": "1",
                "start": 300,
                "end": 400,
                "gene_symbol": "GENE1",
                "exon_id": "E2",
                "strand": "+",
                "transcript_id": "T1",
                "rank": 2,
            },
            # Gene 2, MANE select transcript
            {
                "chromosome": "2",
                "start": 500,
                "end": 550,
                "gene_symbol": "GENE2",
                "exon_id": "E3",
                "strand": "-",
                "transcript_id": "T2",
                "rank": 1,
            },
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "test_exons.bed"
            create_exon_bed_file(
                exons_data, bed_path, transcript_type="canonical", padding=10
            )

            assert bed_path.exists()
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 3

            # Check first exon of GENE1 (0-based, padded)
            fields1 = lines[0].strip().split("\t")
            assert fields1[0] == "1"
            assert fields1[1] == "89"  # 100 - 1 - 10
            assert fields1[2] == "210"  # 200 + 10
            assert "GENE1_T1_exon1" in fields1[3]
            assert fields1[4] == "1000"
            assert fields1[5] == "+"

            # Check second exon of GENE1
            fields2 = lines[1].strip().split("\t")
            assert fields2[0] == "1"
            assert fields2[1] == "289"  # 300 - 1 - 10
            assert fields2[2] == "410"  # 400 + 10
            assert "GENE1_T1_exon2" in fields2[3]

            # Check exon of GENE2 (negative strand)
            fields3 = lines[2].strip().split("\t")
            assert fields3[0] == "2"
            assert fields3[1] == "489"  # 500 - 1 - 10
            assert fields3[2] == "560"  # 550 + 10
            assert "GENE2_T2_exon1" in fields3[3]
            assert fields3[5] == "-"


class TestEnsemblHelpers:
    """Test helper functions within the EnsemblClient."""

    @pytest.fixture
    def ensembl_client(self) -> EnsemblClient:
        """Create an EnsemblClient for testing helper functions."""
        # We don't need network access for these tests
        return EnsemblClient(cache_manager=None)

    def test_calculate_transcript_coverage(self, ensembl_client: EnsemblClient):
        """Test transcript coverage calculation with and without padding."""
        transcript_data = {
            "Exon": [
                {"start": 100, "end": 200},  # length 101
                {"start": 300, "end": 400},  # length 101
            ]
        }
        # Total exon length = 101 + 101 = 202

        # Test without padding
        coverage = ensembl_client.calculate_transcript_coverage(
            transcript_data, padding=0
        )
        assert coverage == 202

        # Test with padding
        coverage_padded = ensembl_client.calculate_transcript_coverage(
            transcript_data, padding=10
        )
        # Expected: 202 + (2 * 10) = 222
        assert coverage_padded == 222

    def test_calculate_transcript_coverage_edge_cases(
        self, ensembl_client: EnsemblClient
    ):
        """Test edge cases for transcript coverage calculation."""
        # No exons - returns 0, not None
        assert ensembl_client.calculate_transcript_coverage({"Exon": []}) == 0
        # Malformed data - returns None
        assert ensembl_client.calculate_transcript_coverage({}) is None
        assert ensembl_client.calculate_transcript_coverage(None) is None

    def test_calculate_transcript_coverage_overlapping_exons(
        self, ensembl_client: EnsemblClient
    ):
        """Test transcript coverage calculation with overlapping exons."""
        transcript_data = {
            "Exon": [
                {"start": 100, "end": 250},  # length 151
                {"start": 200, "end": 300},  # length 101, overlaps with first
            ]
        }
        # Should handle overlapping exons properly
        coverage = ensembl_client.calculate_transcript_coverage(
            transcript_data, padding=0
        )
        # This would depend on the implementation - might be sum of all exons or merged regions
        assert coverage is not None
        assert coverage > 0


class TestSNPDeduplication:
    """Test SNP deduplication functionality."""

    def test_deduplicate_snps_basic(self):
        """Test basic SNP deduplication with duplicate VCF IDs."""
        from custom_panel.engine.pipeline import Pipeline

        # Create a mock pipeline instance to access the deduplication method
        pipeline = Pipeline({"output": {"directory": "/tmp"}})
        _deduplicate_snps = pipeline._deduplicate_snps_in_pipeline

        # Create test data with duplicate VCF IDs from different sources
        df = pd.DataFrame({
            "snp": ["19:11091630:G:T", "19:11091630:G:T", "rs123456"],
            "rsid": [pd.NA, "rs6511720", "rs123456"],
            "source": ["PGS_Catalog", "Manual_SNPs", "Identity_SNPs"],
            "category": ["prs", "manual", "identity"],
            "snp_type": ["prs", "manual_snps", "identity"],
            "hg38_chromosome": ["19", "19", "1"],
            "hg38_start": [11091630, 11091630, 1000000],
            "hg38_end": [11091630, 11091630, 1000000],
        })

        result = _deduplicate_snps(df)

        # Should have 2 unique SNPs now
        assert len(result) == 2

        # Check the deduplicated entry for 19:11091630:G:T
        dup_entry = result[result["snp"] == "19:11091630:G:T"].iloc[0]
        assert dup_entry["rsid"] == "rs6511720"  # Should take first non-null rsID
        assert dup_entry["source"] == "Manual_SNPs; PGS_Catalog"  # Sources merged and sorted
        assert dup_entry["category"] == "manual; prs"  # Categories merged and sorted
        assert dup_entry["snp_type"] == "manual_snps; prs"  # Types merged and sorted
        assert dup_entry["source_count"] == 2  # Two sources

        # Check the non-duplicate entry remains unchanged
        single_entry = result[result["snp"] == "rs123456"].iloc[0]
        assert single_entry["rsid"] == "rs123456"
        assert single_entry["source"] == "Identity_SNPs"
        assert single_entry["source_count"] == 1

    def test_deduplicate_snps_no_duplicates(self):
        """Test deduplication with no duplicate entries."""
        from custom_panel.engine.pipeline import Pipeline

        # Create a mock pipeline instance to access the deduplication method
        pipeline = Pipeline({"output": {"directory": "/tmp"}})
        _deduplicate_snps = pipeline._deduplicate_snps_in_pipeline

        # Create test data with no duplicates
        df = pd.DataFrame({
            "snp": ["19:11091630:G:T", "rs123456", "1:1000000:A:G"],
            "rsid": ["rs6511720", "rs123456", "rs789012"],
            "source": ["PGS_Catalog", "Manual_SNPs", "Identity_SNPs"],
            "category": ["prs", "manual", "identity"],
            "snp_type": ["prs", "manual_snps", "identity"],
        })

        result = _deduplicate_snps(df)

        # Should have same number of SNPs
        assert len(result) == 3

        # All should have source_count = 1
        assert all(result["source_count"] == 1)

    def test_deduplicate_snps_empty_dataframe(self):
        """Test deduplication with empty DataFrame."""
        from custom_panel.engine.pipeline import Pipeline

        # Create a mock pipeline instance to access the deduplication method
        pipeline = Pipeline({"output": {"directory": "/tmp"}})
        _deduplicate_snps = pipeline._deduplicate_snps_in_pipeline

        df = pd.DataFrame()
        result = _deduplicate_snps(df)

        assert result.empty

    def test_deduplicate_snps_missing_snp_column(self):
        """Test deduplication with missing 'snp' column."""
        from custom_panel.engine.pipeline import Pipeline

        # Create a mock pipeline instance to access the deduplication method
        pipeline = Pipeline({"output": {"directory": "/tmp"}})
        _deduplicate_snps = pipeline._deduplicate_snps_in_pipeline

        df = pd.DataFrame({
            "rsid": ["rs123456", "rs789012"],
            "source": ["Source1", "Source2"],
        })

        result = _deduplicate_snps(df)

        # Should return original DataFrame unchanged
        assert len(result) == 2
        assert "snp" not in result.columns

    def test_deduplicate_snps_complex_metadata(self):
        """Test deduplication with complex metadata columns."""
        from custom_panel.engine.pipeline import Pipeline

        # Create a mock pipeline instance to access the deduplication method
        pipeline = Pipeline({"output": {"directory": "/tmp"}})
        _deduplicate_snps = pipeline._deduplicate_snps_in_pipeline

        # Create test data with various metadata columns
        df = pd.DataFrame({
            "snp": ["19:11091630:G:T", "19:11091630:G:T", "19:11091630:G:T"],
            "rsid": [pd.NA, "rs6511720", pd.NA],
            "source": ["PGS_Catalog_Clinical", "Manual_SNPs", "PGS_Catalog_Other"],
            "category": ["prs", "manual", "prs"],
            "snp_type": ["prs", "manual_snps", "prs"],
            "pgs_id": ["PGS000765", pd.NA, "PGS000066"],
            "pgs_name": ["Breast Cancer PRS", pd.NA, "Prostate Cancer PRS"],
            "trait": ["breast_cancer", "pharmacogenomics", "prostate_cancer"],
            "effect_weight": [0.5, pd.NA, 0.3],
            "effect_allele": ["T", "T", "T"],
            "hg38_chromosome": ["19", "19", "19"],
            "hg38_start": [11091630, 11091630, 11091630],
        })

        result = _deduplicate_snps(df)

        # Should have 1 unique SNP
        assert len(result) == 1

        entry = result.iloc[0]
        assert entry["snp"] == "19:11091630:G:T"
        assert entry["rsid"] == "rs6511720"  # First non-null value
        assert entry["source"] == "Manual_SNPs; PGS_Catalog_Clinical; PGS_Catalog_Other"  # Sorted
        assert entry["category"] == "manual; prs"  # Unique categories sorted
        assert entry["pgs_id"] == "PGS000066; PGS000765"  # Merged PGS IDs sorted
        assert entry["pgs_name"] == "Breast Cancer PRS; Prostate Cancer PRS"  # Merged names sorted
        assert entry["trait"] == "breast_cancer; pharmacogenomics; prostate_cancer"  # Sorted
        assert entry["effect_weight"] == 0.5  # First non-null value
        assert entry["effect_allele"] == "T"  # First non-null value
        assert entry["source_count"] == 3

    def test_deduplicate_snps_rsid_priority(self):
        """Test that rsID is properly prioritized (first non-null value)."""
        from custom_panel.engine.pipeline import Pipeline

        # Create a mock pipeline instance to access the deduplication method
        pipeline = Pipeline({"output": {"directory": "/tmp"}})
        _deduplicate_snps = pipeline._deduplicate_snps_in_pipeline

        # Create test data where rsID appears in different positions
        df = pd.DataFrame({
            "snp": ["19:11091630:G:T", "19:11091630:G:T", "19:11091630:G:T"],
            "rsid": [pd.NA, "rs6511720", "rs9999999"],  # Different rsIDs
            "source": ["Source1", "Source2", "Source3"],
            "category": ["cat1", "cat2", "cat3"],
            "snp_type": ["type1", "type2", "type3"],
        })

        result = _deduplicate_snps(df)

        # Should take the first non-null rsID
        assert len(result) == 1
        assert result.iloc[0]["rsid"] == "rs6511720"  # First non-null value

    def test_deduplicate_snps_coordinate_preservation(self):
        """Test that coordinate information is properly preserved."""
        from custom_panel.engine.pipeline import Pipeline

        # Create a mock pipeline instance to access the deduplication method
        pipeline = Pipeline({"output": {"directory": "/tmp"}})
        _deduplicate_snps = pipeline._deduplicate_snps_in_pipeline

        # Create test data with coordinate information
        df = pd.DataFrame({
            "snp": ["19:11091630:G:T", "19:11091630:G:T"],
            "rsid": [pd.NA, "rs6511720"],
            "source": ["Source1", "Source2"],
            "category": ["cat1", "cat2"],
            "snp_type": ["type1", "type2"],
            "hg38_chromosome": ["19", "19"],
            "hg38_start": [11091630, 11091630],
            "hg38_end": [11091630, 11091630],
            "hg38_strand": ["+", "+"],
            "ref_allele": ["G", "G"],
            "alt_allele": ["T", "T"],
        })

        result = _deduplicate_snps(df)

        # Should preserve coordinate information
        assert len(result) == 1
        entry = result.iloc[0]
        assert entry["hg38_chromosome"] == "19"
        assert entry["hg38_start"] == 11091630
        assert entry["hg38_end"] == 11091630
        assert entry["hg38_strand"] == "+"
        assert entry["ref_allele"] == "G"
        assert entry["alt_allele"] == "T"


class TestNewBEDFileGeneration:
    """Test new master BED file generation functions."""

    def test_create_genes_all_bed(self):
        """Test creating BED file with all genes."""
        from custom_panel.core.io import create_genes_all_bed

        # Create test data
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "TP53", "BRCA2"],
                "chromosome": ["chr17", "chr17", "chr13"],
                "gene_start": [43044295, 7668402, 32315086],
                "gene_end": [43125364, 7687538, 32400266],
                "gene_strand": [1, -1, 1],
                "include": [True, False, True],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "genes_all.bed"
            create_genes_all_bed(df, bed_path, padding=0)

            # Should include all genes regardless of inclusion status
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 3  # All genes included

            # Check that file contains all genes
            bed_content = "".join(lines)
            assert "BRCA1" in bed_content
            assert "TP53" in bed_content
            assert "BRCA2" in bed_content

    def test_create_genes_included_bed(self):
        """Test creating BED file with only included genes."""
        from custom_panel.core.io import create_genes_included_bed

        # Create test data
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "TP53", "BRCA2"],
                "chromosome": ["chr17", "chr17", "chr13"],
                "gene_start": [43044295, 7668402, 32315086],
                "gene_end": [43125364, 7687538, 32400266],
                "gene_strand": [1, -1, 1],
                "include": [True, False, True],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "genes_included.bed"
            create_genes_included_bed(df, bed_path, padding=0)

            # Should only include genes with include=True
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 2  # Only BRCA1 and BRCA2 included

            # Check that file contains only included genes
            bed_content = "".join(lines)
            assert "BRCA1" in bed_content
            assert "TP53" not in bed_content
            assert "BRCA2" in bed_content

    def test_create_snps_all_bed(self):
        """Test creating BED file with all SNPs combined."""
        from custom_panel.core.io import create_snps_all_bed

        # Create test SNP data
        snp_data = {
            "identity": pd.DataFrame({
                "snp": ["rs1234", "rs5678"],
                "hg38_chromosome": ["1", "2"],
                "hg38_start": [1000, 2000],
                "hg38_end": [1000, 2000],
            }),
            "ethnicity": pd.DataFrame({
                "snp": ["rs9876", "rs5432"],
                "hg38_chromosome": ["3", "4"],
                "hg38_start": [3000, 4000],
                "hg38_end": [3000, 4000],
            }),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "snps_all.bed"
            create_snps_all_bed(snp_data, bed_path)

            # Should include all SNPs from all categories
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 4  # All SNPs included

            # Check that file contains all SNPs
            bed_content = "".join(lines)
            assert "rs1234" in bed_content
            assert "rs5678" in bed_content
            assert "rs9876" in bed_content
            assert "rs5432" in bed_content

    def test_create_regions_all_bed(self):
        """Test creating BED file with all regions combined."""
        from custom_panel.core.io import create_regions_all_bed

        # Create test regions data
        regions_data = {
            "manual": pd.DataFrame({
                "region_name": ["region1", "region2"],
                "chromosome": ["chr1", "chr2"],
                "start": [1000, 2000],
                "end": [2000, 3000],
            }),
            "stripy": pd.DataFrame({
                "region_name": ["region3", "region4"],
                "chromosome": ["chr3", "chr4"],
                "start": [3000, 4000],
                "end": [4000, 5000],
            }),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "regions_all.bed"
            create_regions_all_bed(regions_data, bed_path)

            # Should include all regions from all categories
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 4  # All regions included

            # Check that file contains all regions
            bed_content = "".join(lines)
            assert "region1" in bed_content
            assert "region2" in bed_content
            assert "region3" in bed_content
            assert "region4" in bed_content

    def test_create_complete_panel_bed(self):
        """Test creating complete panel BED file with only included genes plus SNPs and regions."""
        from custom_panel.core.io import create_complete_panel_bed

        # Create test gene data
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "TP53"],
                "chromosome": ["chr17", "chr17"],
                "gene_start": [43044295, 7668402],
                "gene_end": [43125364, 7687538],
                "gene_strand": [1, -1],
                "include": [True, False],  # Only BRCA1 is included
            }
        )

        # Create test SNP data
        snp_data = {
            "identity": pd.DataFrame({
                "snp": ["rs1234"],
                "hg38_chromosome": ["1"],
                "hg38_start": [1000],
                "hg38_end": [1000],
            }),
        }

        # Create test regions data
        regions_data = {
            "manual": pd.DataFrame({
                "region_name": ["region1"],
                "chromosome": ["chr1"],
                "start": [2000],
                "end": [3000],
            }),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "complete_panel.bed"
            create_complete_panel_bed(df, snp_data, regions_data, bed_path, padding=0)

            # Should include only included genes + SNPs + regions
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 3  # 1 included gene + 1 SNP + 1 region

            # Check that file contains only included genes plus SNPs and regions
            bed_content = "".join(lines)
            assert "BRCA1" in bed_content
            assert "TP53" not in bed_content  # Not included
            assert "rs1234" in bed_content
            assert "region1" in bed_content

    def test_create_complete_panel_bed_with_padding(self):
        """Test creating complete panel BED file with padding."""
        from custom_panel.core.io import create_complete_panel_bed

        # Create test gene data
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1"],
                "chromosome": ["chr17"],
                "gene_start": [43044295],
                "gene_end": [43125364],
                "gene_strand": [1],
                "include": [True],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "complete_panel_padded.bed"
            create_complete_panel_bed(df, None, None, bed_path, padding=100)

            # Check that padding was applied
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 1
            line = lines[0].strip().split("\t")

            # BED format: chr, start, end, name, score, strand, element_type, element_subtype
            assert line[0] == "chr17"
            assert int(line[1]) == 43044295 - 1 - 100  # BED is 0-based with padding
            assert int(line[2]) == 43125364 + 100  # End position with padding
            assert line[3] == "BRCA1"

    def test_empty_data_handling(self):
        """Test handling of empty data for all new BED functions."""
        from custom_panel.core.io import (
            create_complete_panel_bed,
            create_genes_all_bed,
            create_regions_all_bed,
            create_snps_all_bed,
        )

        # Test with empty DataFrame
        empty_df = pd.DataFrame()

        with tempfile.TemporaryDirectory() as tmpdir:
            # Test genes_all_bed with empty data
            bed_path = Path(tmpdir) / "empty_genes_all.bed"
            # Should handle gracefully and not create file or create empty file
            try:
                create_genes_all_bed(empty_df, bed_path, padding=0)
            except ValueError:
                # Expected behavior for missing required columns
                pass

            # Test snps_all_bed with empty data
            bed_path = Path(tmpdir) / "empty_snps_all.bed"
            create_snps_all_bed({}, bed_path)
            # Should handle gracefully

            # Test regions_all_bed with empty data
            bed_path = Path(tmpdir) / "empty_regions_all.bed"
            create_regions_all_bed({}, bed_path)
            # Should handle gracefully

            # Test complete_panel_bed with empty data
            bed_path = Path(tmpdir) / "empty_complete_panel.bed"
            create_complete_panel_bed(empty_df, None, None, bed_path, padding=0)
            # Should handle gracefully

    def test_bed_file_sorting(self):
        """Test that BED files are properly sorted by chromosome and position."""
        from custom_panel.core.io import create_genes_all_bed

        # Create test data with chromosomes in non-natural order
        df = pd.DataFrame(
            {
                "approved_symbol": ["GENE1", "GENE2", "GENE3", "GENE4"],
                "chromosome": ["chr10", "chr2", "chrX", "chr1"],
                "gene_start": [1000, 2000, 3000, 4000],
                "gene_end": [2000, 3000, 4000, 5000],
                "gene_strand": [1, 1, 1, 1],
                "include": [True, True, True, True],
            }
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "sorted_genes.bed"
            create_genes_all_bed(df, bed_path, padding=0)

            # Check that file is sorted naturally (1, 2, 10, X)
            with open(bed_path) as f:
                lines = f.readlines()

            # First line should be chr1 (GENE4)
            assert lines[0].split("\t")[0] == "chr1"
            assert "GENE4" in lines[0]

            # Second line should be chr2 (GENE2)
            assert lines[1].split("\t")[0] == "chr2"
            assert "GENE2" in lines[1]

            # Third line should be chr10 (GENE1)
            assert lines[2].split("\t")[0] == "chr10"
            assert "GENE1" in lines[2]

            # Fourth line should be chrX (GENE3)
            assert lines[3].split("\t")[0] == "chrX"
            assert "GENE3" in lines[3]

    def test_create_complete_panel_exons_bed(self):
        """Test creating complete panel exons BED file with exons from included genes plus SNPs and regions."""
        from custom_panel.core.io import create_complete_panel_exons_bed

        # Create test gene data
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "TP53"],
                "chromosome": ["chr17", "chr17"],
                "gene_start": [43044295, 7668402],
                "gene_end": [43125364, 7687538],
                "gene_strand": [1, -1],
                "include": [True, False],  # Only BRCA1 is included
                "canonical_transcript": ["ENST00000357654", "ENST00000269305"],
            }
        )

        # Create mock transcript data
        transcript_data = {
            "BRCA1": {
                "all_transcripts": [
                    {
                        "id": "ENST00000357654",
                        "Exon": [
                            {"seq_region_name": "17", "start": 43044295, "end": 43044395, "strand": 1, "rank": 1},
                            {"seq_region_name": "17", "start": 43045000, "end": 43045100, "strand": 1, "rank": 2},
                        ]
                    }
                ]
            }
        }

        # Create test SNP data
        snp_data = {
            "identity": pd.DataFrame({
                "snp": ["rs1234"],
                "hg38_chromosome": ["1"],
                "hg38_start": [1000],
                "hg38_end": [1000],
            }),
        }

        # Create test regions data
        regions_data = {
            "manual": pd.DataFrame({
                "region_name": ["region1"],
                "chromosome": ["chr1"],
                "start": [2000],
                "end": [3000],
            }),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "complete_panel_exons.bed"
            create_complete_panel_exons_bed(df, transcript_data, snp_data, regions_data, bed_path, padding=0)

            # Should include exons from included genes + SNPs + regions
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 4  # 2 exons + 1 SNP + 1 region

            # Check that file contains exons from included genes plus SNPs and regions
            bed_content = "".join(lines)
            assert "BRCA1_exon1" in bed_content
            assert "BRCA1_exon2" in bed_content
            assert "rs1234" in bed_content
            assert "region1" in bed_content

    def test_create_complete_panel_genes_bed(self):
        """Test creating complete panel genes BED file with full genomic regions from included genes plus SNPs and regions."""
        from custom_panel.core.io import create_complete_panel_genes_bed

        # Create test gene data
        df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "TP53"],
                "chromosome": ["chr17", "chr17"],
                "gene_start": [43044295, 7668402],
                "gene_end": [43125364, 7687538],
                "gene_strand": [1, -1],
                "include": [True, False],  # Only BRCA1 is included
            }
        )

        # Create test SNP data
        snp_data = {
            "identity": pd.DataFrame({
                "snp": ["rs1234"],
                "hg38_chromosome": ["1"],
                "hg38_start": [1000],
                "hg38_end": [1000],
            }),
        }

        # Create test regions data
        regions_data = {
            "manual": pd.DataFrame({
                "region_name": ["region1"],
                "chromosome": ["chr1"],
                "start": [2000],
                "end": [3000],
            }),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "complete_panel_genes.bed"
            create_complete_panel_genes_bed(df, snp_data, regions_data, bed_path, padding=0)

            # Should include full genomic regions from included genes + SNPs + regions
            with open(bed_path) as f:
                lines = f.readlines()

            assert len(lines) == 3  # 1 included gene + 1 SNP + 1 region

            # Check that file contains full genomic regions from included genes plus SNPs and regions
            bed_content = "".join(lines)
            assert "BRCA1" in bed_content
            assert "TP53" not in bed_content  # Not included
            assert "rs1234" in bed_content
            assert "region1" in bed_content
