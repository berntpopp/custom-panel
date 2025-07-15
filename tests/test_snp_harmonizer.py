"""
Unit tests for SNP harmonizer functionality.
"""

from __future__ import annotations

from unittest.mock import Mock

import pandas as pd
import pytest

from custom_panel.core.snp_harmonizer import SNPHarmonizer


class TestSNPHarmonizer:
    """Test the SNP harmonizer functionality."""

    @pytest.fixture
    def mock_gnomad_client(self):
        """Create a mock gnomAD client."""
        client = Mock()
        client.resolve_rsid.return_value = "rs123456"
        client.liftover_coordinates.return_value = ("1", "1-54321-A-T", 54321, "A", "T")
        client.validate_variant.return_value = True
        return client

    @pytest.fixture
    def mock_ensembl_client(self):
        """Create a mock Ensembl client."""
        client = Mock()
        client.get_variations_batch.return_value = {}
        client.extract_coordinates_from_variation.return_value = None
        return client

    @pytest.fixture
    def harmonizer(self, mock_gnomad_client, mock_ensembl_client):
        """Create a test SNP harmonizer."""
        config = {
            "prefer_rsid": True,
            "fallback_to_coordinates": True,
        }
        return SNPHarmonizer(mock_gnomad_client, mock_ensembl_client, config)

    def test_harmonizer_initialization(
        self, harmonizer, mock_gnomad_client, mock_ensembl_client
    ):
        """Test harmonizer initialization."""
        assert harmonizer.gnomad_client == mock_gnomad_client
        assert harmonizer.ensembl_client == mock_ensembl_client
        assert harmonizer.prefer_rsid is True
        assert harmonizer.fallback_to_coordinates is True

    def test_harmonize_empty_dataframe(self, harmonizer):
        """Test harmonizing empty DataFrame."""
        empty_df = pd.DataFrame()
        result = harmonizer.harmonize_snp_batch(empty_df)
        assert result.empty

    def test_harmonize_snp_with_rsid(self, harmonizer, mock_gnomad_client):
        """Test harmonizing SNP that already has rsID."""
        snp_data = pd.DataFrame(
            [
                {
                    "snp": "rs123456",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "A/T",
                }
            ]
        )

        # Mock liftover to return hg19 coordinates
        mock_gnomad_client.liftover_coordinates.return_value = (
            "1",
            "1-54321-A-T",
            54321,
            "A",
            "T",
        )
        mock_gnomad_client.validate_variant.return_value = True

        result = harmonizer.harmonize_snp_batch(snp_data)

        assert len(result) == 1
        assert result.iloc[0]["snp"] == "rs123456"
        assert result.iloc[0]["rsid"] == "rs123456"
        assert result.iloc[0]["harmonized"]

    def test_harmonize_snp_resolve_rsid(self, harmonizer, mock_gnomad_client):
        """Test harmonizing SNP by resolving rsID from coordinates."""
        snp_data = pd.DataFrame(
            [
                {
                    "snp": "coordinate_id",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "A/T",
                }
            ]
        )

        mock_gnomad_client.resolve_rsid.return_value = "rs789012"

        result = harmonizer.harmonize_snp_batch(snp_data)

        assert len(result) == 1
        assert result.iloc[0]["snp"] == "rs789012"
        assert result.iloc[0]["rsid"] == "rs789012"
        assert result.iloc[0]["original_id"] == "coordinate_id"

    def test_harmonize_snp_fallback_to_coordinates(
        self, harmonizer, mock_gnomad_client
    ):
        """Test harmonizing SNP falling back to coordinate-based ID."""
        snp_data = pd.DataFrame(
            [
                {
                    "snp": "unknown_id",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "A/T",
                }
            ]
        )

        mock_gnomad_client.resolve_rsid.return_value = None
        mock_gnomad_client.liftover_coordinates.return_value = (
            "1",
            "1-54321-A-T",
            54321,
            "A",
            "T",
        )
        mock_gnomad_client.validate_variant.return_value = True

        result = harmonizer.harmonize_snp_batch(snp_data)

        assert len(result) == 1
        assert result.iloc[0]["snp"] == "1:12345:A:T"
        assert result.iloc[0]["coordinate_based_id"]
        assert result.iloc[0]["original_id"] == "unknown_id"

    def test_harmonize_snp_liftover_coordinates(self, harmonizer, mock_gnomad_client):
        """Test harmonizing SNP with coordinate liftover."""
        snp_data = pd.DataFrame(
            [
                {
                    "snp": "rs123456",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "A/T",
                }
            ]
        )

        mock_gnomad_client.liftover_coordinates.return_value = (
            "1",
            "1-54321-A-T",
            54321,
            "A",
            "T",
        )

        result = harmonizer.harmonize_snp_batch(snp_data)

        assert len(result) == 1
        assert result.iloc[0]["hg19_chromosome"] == "1"
        assert result.iloc[0]["hg19_start"] == 54321
        assert result.iloc[0]["hg19_allele_string"] == "A/T"

    # NOTE: Validation tests are disabled as gnomAD validation is no longer used for rsID processing
    # def test_harmonize_snp_validation_success(self, harmonizer, mock_gnomad_client):
    #     """Test harmonizing SNP with successful validation."""
    #     # This test is disabled as validation is no longer performed for rsID processing
    #     pass

    # def test_harmonize_snp_validation_failure(self, harmonizer, mock_gnomad_client):
    #     """Test harmonizing SNP with validation failure."""
    #     # This test is disabled as validation is no longer performed for rsID processing
    #     pass

    def test_harmonize_snp_error_handling(self, harmonizer, mock_gnomad_client):
        """Test error handling during SNP harmonization."""
        snp_data = pd.DataFrame(
            [
                {
                    "snp": "coordinate_id",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "A/T",
                }
            ]
        )

        mock_gnomad_client.resolve_rsid.side_effect = Exception("API Error")

        result = harmonizer.harmonize_snp_batch(snp_data)

        assert len(result) == 1
        assert "harmonization_error" in result.iloc[0]
        assert result.iloc[0]["harmonization_error"] == "API Error"

    def test_has_coordinates_methods(self, harmonizer):
        """Test coordinate detection methods."""
        # Test with hg38 coordinates
        snp_info_hg38 = {"hg38_chromosome": "1", "hg38_start": 12345}
        assert harmonizer._has_coordinates(snp_info_hg38) is True
        assert harmonizer._has_build_coordinates(snp_info_hg38, "hg38") is True
        assert harmonizer._has_build_coordinates(snp_info_hg38, "hg19") is False

        # Test with hg19 coordinates
        snp_info_hg19 = {"hg19_chromosome": "1", "hg19_start": 12345}
        assert harmonizer._has_coordinates(snp_info_hg19) is True
        assert harmonizer._has_build_coordinates(snp_info_hg19, "hg19") is True
        assert harmonizer._has_build_coordinates(snp_info_hg19, "hg38") is False

        # Test with generic coordinates
        snp_info_generic = {"chromosome": "1", "position": 12345}
        assert harmonizer._has_coordinates(snp_info_generic) is True

        # Test without coordinates
        snp_info_no_coords = {"snp": "rs123456"}
        assert harmonizer._has_coordinates(snp_info_no_coords) is False

    def test_extract_alleles(self, harmonizer):
        """Test allele extraction methods."""
        # Test with allele string
        snp_info = {"hg38_allele_string": "A/T"}
        ref_allele = harmonizer._extract_ref_allele(snp_info, "hg38")
        alt_allele = harmonizer._extract_alt_allele(snp_info, "hg38")
        assert ref_allele == "A"
        assert alt_allele == "T"

        # Test with separate allele columns
        snp_info = {"ref_allele": "G", "alt_allele": "C"}
        ref_allele = harmonizer._extract_ref_allele(snp_info, "hg38")
        alt_allele = harmonizer._extract_alt_allele(snp_info, "hg38")
        assert ref_allele == "G"
        assert alt_allele == "C"

    def test_create_coordinate_id(self, harmonizer):
        """Test coordinate ID creation."""
        # Test with hg38 coordinates
        snp_info_hg38 = {
            "hg38_chromosome": "1",
            "hg38_start": 12345,
            "hg38_allele_string": "A/T",
        }
        coord_id = harmonizer._create_coordinate_id(snp_info_hg38)
        assert coord_id == "1:12345:A:T"

        # Test with hg19 coordinates
        snp_info_hg19 = {
            "hg19_chromosome": "1",
            "hg19_start": 12345,
            "hg19_allele_string": "A/T",
        }
        coord_id = harmonizer._create_coordinate_id(snp_info_hg19)
        assert coord_id == "1:12345:A:T"

        # Test without coordinates
        snp_info_no_coords = {"snp": "rs123456"}
        coord_id = harmonizer._create_coordinate_id(snp_info_no_coords)
        assert coord_id is None

    def test_extract_coordinates_from_id(self, harmonizer):
        """Test coordinate extraction from ID."""
        snp_info = {"snp": "1:12345:A:T"}

        harmonizer._extract_coordinates_from_id(snp_info)

        assert snp_info["hg38_chromosome"] == "1"
        assert snp_info["hg38_start"] == 12345
        assert snp_info["hg38_end"] == 12345
        assert snp_info["hg38_strand"] == 1
        assert snp_info["hg38_allele_string"] == "A/T"

    def test_harmonizer_stats(self, harmonizer):
        """Test harmonizer statistics."""
        # Initial stats
        stats = harmonizer.get_stats()
        assert stats["total_processed"] == 0
        assert stats["rsid_resolved"] == 0
        assert stats["errors"] == 0

        # Reset stats
        harmonizer.reset_stats()
        stats = harmonizer.get_stats()
        assert stats["total_processed"] == 0

    def test_harmonizer_config_defaults(self, mock_gnomad_client, mock_ensembl_client):
        """Test harmonizer with default configuration."""
        harmonizer = SNPHarmonizer(mock_gnomad_client, mock_ensembl_client)

        assert harmonizer.prefer_rsid is True
        assert harmonizer.fallback_to_coordinates is True

    def test_harmonizer_custom_config(self, mock_gnomad_client, mock_ensembl_client):
        """Test harmonizer with custom configuration."""
        config = {
            "prefer_rsid": False,
            "fallback_to_coordinates": False,
        }
        harmonizer = SNPHarmonizer(mock_gnomad_client, mock_ensembl_client, config)

        assert harmonizer.prefer_rsid is False
        assert harmonizer.fallback_to_coordinates is False

    def test_harmonize_batch_with_mixed_data(self, harmonizer, mock_gnomad_client):
        """Test harmonizing batch with mixed data types."""
        snp_data = pd.DataFrame(
            [
                {
                    "snp": "rs123456",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "A/T",
                },
                {
                    "snp": "coordinate_id",
                    "hg38_chromosome": "2",
                    "hg38_start": 67890,
                    "hg38_allele_string": "G/C",
                },
                {
                    "snp": "1:54321:T:G",
                    "hg19_chromosome": "1",
                    "hg19_start": 54321,
                    "hg19_allele_string": "T/G",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "T/G",
                },
            ]
        )

        mock_gnomad_client.resolve_rsid.return_value = "rs789012"
        mock_gnomad_client.liftover_coordinates.return_value = (
            "1",
            "1-12345-T-G",
            12345,
            "T",
            "G",
        )
        mock_gnomad_client.validate_variant.return_value = True

        result = harmonizer.harmonize_snp_batch(snp_data)

        assert len(result) == 3
        assert all(row["harmonized"] for _, row in result.iterrows())
        assert harmonizer.get_stats()["total_processed"] == 3
