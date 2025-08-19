"""
Tests for PharmGKB pharmacogenomics SNPs source fetcher.
"""

import pickle
import tempfile
import zipfile
from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import Mock, patch, mock_open

import pandas as pd
import pytest
import requests

from custom_panel.sources_snp.pharmacogenomics_snps import (
    _download_and_process_pharmgkb_data,
    _download_pharmgkb_variants,
    _extract_rsid,
    _filter_clinically_relevant_variants,
    _is_cache_valid,
    _transform_to_snp_format,
    fetch_pharmacogenomics_snps,
    get_pharmacogenomics_snps_summary,
)


class TestFetchPharmacogenomicsSnps:
    """Test main fetch function."""

    def test_snp_processing_disabled(self):
        """Test when SNP processing is disabled."""
        config = {"snp_processing": {"enabled": False}}
        result = fetch_pharmacogenomics_snps(config)
        assert result is None

    def test_pharmacogenomics_disabled(self):
        """Test when pharmacogenomics SNPs are disabled."""
        config = {
            "snp_processing": {
                "enabled": True,
                "pharmacogenomics": {"enabled": False},
            }
        }
        result = fetch_pharmacogenomics_snps(config)
        assert result is None

    @patch("custom_panel.sources_snp.pharmacogenomics_snps._download_and_process_pharmgkb_data")
    @patch("custom_panel.sources_snp.pharmacogenomics_snps._filter_clinically_relevant_variants")
    @patch("custom_panel.sources_snp.pharmacogenomics_snps._transform_to_snp_format")
    def test_successful_fetch(self, mock_transform, mock_filter, mock_download):
        """Test successful SNP fetching."""
        # Setup mocks
        mock_variants_df = pd.DataFrame({
            "Variant ID": ["PA166156302"],
            "Variant Name": ["rs4244285"],
            "Gene Symbols": ["CYP2C19"],
            "Guideline Annotation count": [2],
            "Level 1/2 Clinical Annotation count": [1],
        })
        mock_download.return_value = mock_variants_df

        mock_filtered_df = mock_variants_df.copy()
        mock_filter.return_value = mock_filtered_df

        mock_snp_df = pd.DataFrame({
            "snp": ["rs4244285"],
            "rsid": ["rs4244285"],
            "source": ["PharmGKB"],
            "category": ["pharmacogenomics"],
            "gene": ["CYP2C19"],
        })
        mock_transform.return_value = mock_snp_df

        config = {
            "snp_processing": {
                "enabled": True,
                "pharmacogenomics": {"enabled": True},
            }
        }

        result = fetch_pharmacogenomics_snps(config)

        assert result is not None
        assert len(result) == 1
        assert result.iloc[0]["rsid"] == "rs4244285"
        mock_download.assert_called_once()
        mock_filter.assert_called_once()
        mock_transform.assert_called_once()

    @patch("custom_panel.sources_snp.pharmacogenomics_snps._download_and_process_pharmgkb_data")
    def test_download_failure(self, mock_download):
        """Test handling of download failure."""
        mock_download.return_value = None

        config = {
            "snp_processing": {
                "enabled": True,
                "pharmacogenomics": {"enabled": True},
            }
        }

        result = fetch_pharmacogenomics_snps(config)
        assert result is None

    @patch("custom_panel.sources_snp.pharmacogenomics_snps._download_and_process_pharmgkb_data")
    def test_exception_handling(self, mock_download):
        """Test exception handling in main fetch function."""
        mock_download.side_effect = Exception("Test error")

        config = {
            "snp_processing": {
                "enabled": True,
                "pharmacogenomics": {"enabled": True},
            }
        }

        result = fetch_pharmacogenomics_snps(config)
        assert result is None


class TestDownloadAndProcessPharmgkbData:
    """Test download and caching functionality."""

    def test_cache_hit(self):
        """Test cache hit scenario."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir)
            processed_cache = cache_dir / "processed_variants.pkl"

            # Create mock cached data
            mock_df = pd.DataFrame({"test": ["data"]})
            with open(processed_cache, "wb") as f:
                pickle.dump(mock_df, f)

            config = {
                "cache_dir": str(cache_dir),
                "cache_ttl_days": 7,
            }

            with patch("custom_panel.sources_snp.pharmacogenomics_snps._is_cache_valid", return_value=True):
                result = _download_and_process_pharmgkb_data(config)

            assert result is not None
            assert result.equals(mock_df)

    @patch("custom_panel.sources_snp.pharmacogenomics_snps._download_pharmgkb_variants")
    def test_cache_miss(self, mock_download):
        """Test cache miss scenario."""
        mock_df = pd.DataFrame({"test": ["data"]})
        mock_download.return_value = mock_df

        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir)
            config = {
                "cache_dir": str(cache_dir),
                "cache_ttl_days": 7,
                "download_url": "https://test.com/variants.zip",
            }

            with patch("custom_panel.sources_snp.pharmacogenomics_snps._is_cache_valid", return_value=False):
                result = _download_and_process_pharmgkb_data(config)

            assert result is not None
            assert result.equals(mock_df)
            mock_download.assert_called_once()

            # Check that data was cached
            processed_cache = cache_dir / "processed_variants.pkl"
            assert processed_cache.exists()


class TestDownloadPharmgkbVariants:
    """Test PharmGKB data download functionality."""

    @patch("custom_panel.sources_snp.pharmacogenomics_snps.requests.get")
    @patch("custom_panel.sources_snp.pharmacogenomics_snps.zipfile.ZipFile")
    @patch("pandas.read_csv")
    def test_successful_download(self, mock_read_csv, mock_zipfile, mock_get):
        """Test successful download and extraction."""
        # Mock response
        mock_response = Mock()
        mock_response.headers = {"content-length": "1000"}
        mock_response.iter_content.return_value = [b"test data"]
        mock_get.return_value = mock_response

        # Mock ZIP file
        mock_zip = Mock()
        mock_zip.namelist.return_value = ["variants.tsv"]
        mock_zip.open.return_value.__enter__ = Mock(return_value=Mock())
        mock_zip.open.return_value.__exit__ = Mock(return_value=None)
        mock_zipfile.return_value.__enter__ = Mock(return_value=mock_zip)
        mock_zipfile.return_value.__exit__ = Mock(return_value=None)

        # Mock pandas read_csv
        mock_df = pd.DataFrame({
            "Variant ID": ["PA166156302"],
            "Variant Name": ["rs4244285"],
        })
        mock_read_csv.return_value = mock_df

        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir)
            config = {"timeout": 300}

            result = _download_pharmgkb_variants(
                "https://test.com/variants.zip", cache_dir, config
            )

            assert result is not None
            assert len(result) == 1
            mock_get.assert_called_once()
            mock_response.raise_for_status.assert_called_once()

    @patch("custom_panel.sources_snp.pharmacogenomics_snps.requests.get")
    def test_download_request_failure(self, mock_get):
        """Test request failure handling."""
        mock_get.side_effect = requests.RequestException("Network error")

        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir)
            config = {"timeout": 300}

            result = _download_pharmgkb_variants(
                "https://test.com/variants.zip", cache_dir, config
            )

            assert result is None

    @patch("custom_panel.sources_snp.pharmacogenomics_snps.requests.get")
    def test_zip_file_corruption(self, mock_get):
        """Test handling of corrupted ZIP files."""
        # Mock successful download but bad ZIP
        mock_response = Mock()
        mock_response.headers = {"content-length": "1000"}
        mock_response.iter_content.return_value = [b"corrupted data"]
        mock_get.return_value = mock_response

        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir)
            config = {"timeout": 300}

            with patch("custom_panel.sources_snp.pharmacogenomics_snps.zipfile.ZipFile", 
                      side_effect=zipfile.BadZipFile("Corrupted ZIP")):
                result = _download_pharmgkb_variants(
                    "https://test.com/variants.zip", cache_dir, config
                )

            assert result is None


class TestFilterClinicallyRelevantVariants:
    """Test clinical relevance filtering."""

    def test_or_filtering(self):
        """Test OR logic filtering."""
        df = pd.DataFrame({
            "Guideline Annotation count": [0, 1, 0, 2],
            "Level 1/2 Clinical Annotation count": [1, 0, 0, 1],
            "Variant ID": ["PA1", "PA2", "PA3", "PA4"],
        })

        config = {
            "filters": {
                "min_guideline_annotations": 1,
                "min_level12_clinical_annotations": 1,
                "filter_logic": "OR",
            }
        }

        result = _filter_clinically_relevant_variants(df, config)

        # Should include PA1 (level1/2 >= 1), PA2 (guideline >= 1), PA4 (both)
        # Should exclude PA3 (neither condition met)
        assert len(result) == 3
        assert "PA3" not in result["Variant ID"].values

    def test_and_filtering(self):
        """Test AND logic filtering."""
        df = pd.DataFrame({
            "Guideline Annotation count": [0, 1, 0, 2],
            "Level 1/2 Clinical Annotation count": [1, 0, 0, 1],
            "Variant ID": ["PA1", "PA2", "PA3", "PA4"],
        })

        config = {
            "filters": {
                "min_guideline_annotations": 1,
                "min_level12_clinical_annotations": 1,
                "filter_logic": "AND",
            }
        }

        result = _filter_clinically_relevant_variants(df, config)

        # Only PA4 meets both conditions
        assert len(result) == 1
        assert result.iloc[0]["Variant ID"] == "PA4"

    def test_missing_columns(self):
        """Test handling of missing required columns."""
        df = pd.DataFrame({
            "Variant ID": ["PA1", "PA2"],
            "Other Column": [1, 2],
        })

        config = {"filters": {}}

        result = _filter_clinically_relevant_variants(df, config)
        assert len(result) == 0  # Should return empty DataFrame

    def test_non_numeric_values(self):
        """Test handling of non-numeric annotation counts."""
        df = pd.DataFrame({
            "Guideline Annotation count": ["N/A", "1", "not_a_number"],
            "Level 1/2 Clinical Annotation count": [1, "invalid", 2],
            "Variant ID": ["PA1", "PA2", "PA3"],
        })

        config = {
            "filters": {
                "min_guideline_annotations": 1,
                "min_level12_clinical_annotations": 1,
                "filter_logic": "OR",
            }
        }

        result = _filter_clinically_relevant_variants(df, config)

        # PA1: guideline=0 (converted from N/A), level1/2=1 -> included
        # PA2: guideline=1, level1/2=0 (converted from invalid) -> included  
        # PA3: guideline=0 (converted), level1/2=2 -> included
        assert len(result) == 3


class TestTransformToSnpFormat:
    """Test transformation to standard SNP format."""

    def test_successful_transformation(self):
        """Test successful transformation."""
        df = pd.DataFrame({
            "Variant ID": ["PA166156302", "PA166156303"],
            "Variant Name": ["rs4244285", "rs1234567"],
            "Gene Symbols": ["CYP2C19", "CYP2D6"],
            "Location": ["chr10:94842866", "chr22:42522500"],
            "Clinical Annotation count": [5, 3],
            "Guideline Annotation count": [2, 1],
            "Level 1/2 Clinical Annotation count": [3, 2],
        })

        result = _transform_to_snp_format(df)

        assert len(result) == 2
        assert all(result["category"] == "pharmacogenomics")
        assert all(result["source"] == "PharmGKB")
        assert "rs4244285" in result["rsid"].values
        assert "rs1234567" in result["rsid"].values

    def test_duplicate_handling(self):
        """Test handling of duplicate rsIDs."""
        df = pd.DataFrame({
            "Variant ID": ["PA1", "PA2"],
            "Variant Name": ["rs4244285", "rs4244285"],  # Same rsID
            "Gene Symbols": ["CYP2C19", "CYP2C19"],
            "Location": ["chr10:94842866", "chr10:94842866"],
            "Clinical Annotation count": [3, 5],  # Different annotation counts
            "Guideline Annotation count": [1, 2],
            "Level 1/2 Clinical Annotation count": [2, 3],
        })

        result = _transform_to_snp_format(df)

        # Should keep only one, the one with higher clinical annotation count
        assert len(result) == 1
        assert result.iloc[0]["clinical_annotation_count"] == 5
        assert result.iloc[0]["pharmgkb_id"] == "PA2"

    def test_invalid_variant_names(self):
        """Test handling of variants without valid rsIDs."""
        df = pd.DataFrame({
            "Variant ID": ["PA1", "PA2", "PA3"],
            "Variant Name": ["rs4244285", "invalid_name", None],
            "Gene Symbols": ["CYP2C19", "CYP2D6", "TPMT"],
            "Location": ["chr10:94842866", "chr22:42522500", "chr6:18130918"],
            "Clinical Annotation count": [5, 3, 2],
            "Guideline Annotation count": [2, 1, 1],
            "Level 1/2 Clinical Annotation count": [3, 2, 1],
        })

        result = _transform_to_snp_format(df)

        # Only the first variant should be included
        assert len(result) == 1
        assert result.iloc[0]["rsid"] == "rs4244285"


class TestExtractRsid:
    """Test rsID extraction from Variant Name."""

    def test_simple_rsid(self):
        """Test extraction of simple rsID."""
        assert _extract_rsid("rs4244285") == "rs4244285"

    def test_multiple_rsids(self):
        """Test extraction from multiple rsIDs."""
        assert _extract_rsid("rs4244285; rs1234567") == "rs4244285"  # First one

    def test_mixed_content(self):
        """Test extraction from mixed content."""
        assert _extract_rsid("some text rs4244285 more text") == "rs4244285"

    def test_case_insensitive(self):
        """Test case insensitive extraction."""
        assert _extract_rsid("RS4244285") == "rs4244285"

    def test_no_rsid(self):
        """Test when no rsID is present."""
        assert _extract_rsid("no rsid here") is None
        assert _extract_rsid("1234567") is None
        assert _extract_rsid("") is None
        assert _extract_rsid(None) is None

    def test_invalid_input(self):
        """Test handling of invalid inputs."""
        assert _extract_rsid(123) is None  # Non-string input


class TestIsCacheValid:
    """Test cache validation."""

    def test_cache_valid(self):
        """Test valid cache file."""
        with tempfile.NamedTemporaryFile() as temp_file:
            cache_path = Path(temp_file.name)
            
            # File exists and is recent
            assert _is_cache_valid(cache_path, ttl_days=7) is True

    def test_cache_expired(self):
        """Test expired cache file."""
        with tempfile.NamedTemporaryFile() as temp_file:
            cache_path = Path(temp_file.name)
            
            # Mock old file
            old_time = datetime.now() - timedelta(days=10)
            with patch("custom_panel.sources_snp.pharmacogenomics_snps.datetime") as mock_dt:
                mock_dt.fromtimestamp.return_value = old_time
                mock_dt.now.return_value = datetime.now()
                
                assert _is_cache_valid(cache_path, ttl_days=7) is False

    def test_cache_nonexistent(self):
        """Test non-existent cache file."""
        cache_path = Path("/nonexistent/file.pkl")
        assert _is_cache_valid(cache_path, ttl_days=7) is False


class TestGetPharmacogenomicsSnpsSummary:
    """Test summary statistics generation."""

    def test_empty_dataframe(self):
        """Test summary for empty DataFrame."""
        df = pd.DataFrame()
        summary = get_pharmacogenomics_snps_summary(df)
        assert summary["total_snps"] == 0

    def test_complete_summary(self):
        """Test complete summary generation."""
        df = pd.DataFrame({
            "rsid": ["rs4244285", "rs1234567", "rs4244285"],  # One duplicate
            "source": ["PharmGKB", "PharmGKB", "PharmGKB"],
            "category": ["pharmacogenomics", "pharmacogenomics", "pharmacogenomics"],
            "gene": ["CYP2C19", "CYP2D6;TPMT", "CYP2C19"],
            "clinical_annotation_count": [5, 3, 4],
            "guideline_annotation_count": [2, 0, 1],
            "level12_clinical_annotation_count": [3, 2, 2],
        })

        summary = get_pharmacogenomics_snps_summary(df)

        assert summary["total_snps"] == 3
        assert summary["unique_rsids"] == 2
        assert summary["unique_genes"] == 3  # CYP2C19, CYP2D6, TPMT
        assert summary["with_guidelines"] == 2  # Two have guideline_annotation_count > 0
        assert summary["with_level12_clinical"] == 3  # All have level12_clinical_annotation_count > 0

        # Check statistics
        assert summary["clinical_annotation_count_stats"]["mean"] == 4.0
        assert summary["clinical_annotation_count_stats"]["max"] == 5
        assert summary["clinical_annotation_count_stats"]["min"] == 3

    def test_gene_parsing(self):
        """Test parsing of multiple genes."""
        df = pd.DataFrame({
            "rsid": ["rs1", "rs2", "rs3"],
            "source": ["PharmGKB", "PharmGKB", "PharmGKB"],
            "category": ["pharmacogenomics", "pharmacogenomics", "pharmacogenomics"],
            "gene": ["CYP2C19", "CYP2D6;TPMT", "CYP2C19,CYP3A4|UGT1A1"],
        })

        summary = get_pharmacogenomics_snps_summary(df)

        # Should identify CYP2C19 (appears twice), CYP2D6, TPMT, CYP3A4, UGT1A1
        assert summary["unique_genes"] == 5
        assert "CYP2C19" in summary["gene_breakdown"]
        assert summary["gene_breakdown"]["CYP2C19"] == 2  # Appears in rs1 and rs3


# Integration test for the whole pipeline
class TestPharmacogenomicsIntegration:
    """Integration tests for the complete pharmacogenomics pipeline."""

    @patch("custom_panel.sources_snp.pharmacogenomics_snps.requests.get")
    @patch("custom_panel.sources_snp.pharmacogenomics_snps.zipfile.ZipFile")
    @patch("pandas.read_csv")
    def test_end_to_end_pipeline(self, mock_read_csv, mock_zipfile, mock_get):
        """Test the complete pipeline from download to SNP format."""
        # Mock network response
        mock_response = Mock()
        mock_response.headers = {"content-length": "1000"}
        mock_response.iter_content.return_value = [b"test data"]
        mock_get.return_value = mock_response

        # Mock ZIP file
        mock_zip = Mock()
        mock_zip.namelist.return_value = ["variants.tsv"]
        mock_zip.open.return_value.__enter__ = Mock(return_value=Mock())
        mock_zip.open.return_value.__exit__ = Mock(return_value=None)
        mock_zipfile.return_value.__enter__ = Mock(return_value=mock_zip)
        mock_zipfile.return_value.__exit__ = Mock(return_value=None)

        # Mock PharmGKB data
        mock_pharmgkb_data = pd.DataFrame({
            "Variant ID": ["PA166156302", "PA166156303", "PA166156304"],
            "Variant Name": ["rs4244285", "rs1234567", "no_rsid_variant"],
            "Gene IDs": ["PA128", "PA129", "PA130"],
            "Gene Symbols": ["CYP2C19", "CYP2D6", "TPMT"],
            "Location": ["NC_000010.11:94842866", "NC_000022.11:42522500", "NC_000006.12:18130918"],
            "Variant Annotation count": [10, 8, 5],
            "Clinical Annotation count": [5, 3, 1],
            "Level 1/2 Clinical Annotation count": [3, 2, 0],
            "Guideline Annotation count": [2, 0, 1],
            "Label Annotation count": [1, 1, 0],
            "Synonyms": ["synonym1", "synonym2", "synonym3"],
        })
        mock_read_csv.return_value = mock_pharmgkb_data

        config = {
            "snp_processing": {
                "enabled": True,
                "pharmacogenomics": {
                    "enabled": True,
                    "cache_dir": "/tmp/test_cache",
                    "cache_ttl_days": 7,
                    "download_url": "https://api.pharmgkb.org/v1/download/file/data/variants.zip",
                    "timeout": 300,
                    "filters": {
                        "min_guideline_annotations": 1,
                        "min_level12_clinical_annotations": 1,
                        "filter_logic": "OR",
                    },
                },
            }
        }

        with patch("custom_panel.sources_snp.pharmacogenomics_snps._is_cache_valid", return_value=False):
            result = fetch_pharmacogenomics_snps(config)

        # Should get 2 SNPs (rs4244285 and rs1234567)
        # PA166156304 should be filtered out (no valid rsID)
        # PA166156303 should be included despite guideline_count=0 because level12_clinical_count=2 >= 1
        assert result is not None
        assert len(result) == 2
        
        rsids = set(result["rsid"])
        assert "rs4244285" in rsids
        assert "rs1234567" in rsids
        
        # Check data integrity
        assert all(result["category"] == "pharmacogenomics")
        assert all(result["source"] == "PharmGKB")
        
        # Check that filtering worked correctly
        cyp2c19_row = result[result["gene"] == "CYP2C19"].iloc[0]
        assert cyp2c19_row["guideline_annotation_count"] == 2
        assert cyp2c19_row["level12_clinical_annotation_count"] == 3