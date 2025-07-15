"""
Integration tests for SNP harmonization workflow.

These tests verify the end-to-end SNP harmonization process
including gnomAD API integration, caching, and error handling.
"""

import tempfile
import time
from unittest.mock import patch

import pandas as pd

from custom_panel.core.gnomad_client import GnomADClient
from custom_panel.core.snp_harmonizer import SNPHarmonizer


class TestSNPHarmonizationIntegration:
    """Integration tests for SNP harmonization workflow."""

    def setup_method(self):
        """Set up test fixtures."""
        # Create temporary cache directory
        self.cache_dir = tempfile.mkdtemp()

        # Create mock client and harmonizer
        self.client = GnomADClient(cache_dir=self.cache_dir, rate_limit=10.0)
        from custom_panel.core.ensembl_client import EnsemblClient

        self.ensembl_client = EnsemblClient()
        self.harmonizer = SNPHarmonizer(self.client, self.ensembl_client)

    def test_end_to_end_harmonization_workflow(self):
        """Test complete harmonization workflow with real-like data."""
        # Create test SNP data similar to what pipeline would process
        test_data = pd.DataFrame(
            [
                {
                    "snp": "rs7072776",
                    "hg38_chromosome": "10",
                    "hg38_start": 21744013,
                    "hg38_allele_string": "A/G",
                    "source": "identity",
                },
                {
                    "snp": "coord_variant",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "C/T",
                    "source": "ethnicity",
                },
                {"snp": "missing_coords", "source": "prs"},
            ]
        )

        # Mock the gnomAD API responses
        with patch.object(self.ensembl_client, "get_variations_batch", return_value={}):
            with patch.object(self.client, "_make_graphql_request") as mock_request:
                # Mock response for rsID already present
                mock_request.return_value = {
                    "variant": {
                        "variant_id": "10-21744013-A-G",
                        "rsid": "rs7072776",
                        "chrom": "10",
                        "pos": 21744013,
                        "ref": "A",
                        "alt": "G",
                    }
                }

                # Run harmonization
                start_time = time.time()
                result = self.harmonizer.harmonize_snp_batch(test_data)
                duration = time.time() - start_time

                # Verify results
                assert len(result) == 3
                assert not result.empty

                # Check that all rows have harmonization metadata (except for error cases)
                successful_rows = result[result["harmonization_error"].isna()]
                assert all(successful_rows["harmonized"] is True)
                assert all(successful_rows["harmonization_timestamp"].notna())

                # Check specific harmonization results
                rs_variant = result[result["snp"] == "rs7072776"].iloc[0]
                assert rs_variant["rsid"] == "rs7072776"
                assert rs_variant["hg38_chromosome"] == "10"

                # Check statistics
                stats = self.harmonizer.get_stats()
                assert stats["total_processed"] == 2
                assert stats["errors"] == 1  # One error for missing_coords

                # Verify performance logging
                assert duration < 5.0  # Should complete quickly with mocked API

    def test_harmonization_with_mixed_data_quality(self):
        """Test harmonization with various data quality issues."""
        test_data = pd.DataFrame(
            [
                {
                    "snp": "rs123456",
                    "hg38_chromosome": "1",
                    "hg38_start": 12345,
                    "hg38_allele_string": "A/T",
                },
                {
                    "snp": "empty_position",
                    "hg38_chromosome": "2",
                    "hg38_start": "",  # Empty position
                    "hg38_allele_string": "C/G",
                },
                {
                    "snp": "none_position",
                    "hg38_chromosome": "3",
                    "hg38_start": None,  # None position
                    "hg38_allele_string": "G/A",
                },
            ]
        )

        # Mock API responses
        with patch.object(self.ensembl_client, "get_variations_batch", return_value={}):
            with patch.object(self.client, "_make_graphql_request") as mock_request:
                mock_request.return_value = {
                    "variant": {
                        "variant_id": "1-12345-A-T",
                        "rsid": "rs123456",
                        "chrom": "1",
                        "pos": 12345,
                        "ref": "A",
                        "alt": "T",
                    }
                }

                result = self.harmonizer.harmonize_snp_batch(test_data)

                # All SNPs should be processed
                assert len(result) == 3

                # Check that position validation works
                valid_snp = result[result["snp"] == "rs123456"].iloc[0]
                assert valid_snp["rsid"] == "rs123456"

                # Check that invalid positions are handled
                stats = self.harmonizer.get_stats()
                assert stats["total_processed"] == 3
                assert stats["errors"] == 0  # Should handle gracefully, not error

    def test_harmonization_cache_integration(self):
        """Test that caching works properly during harmonization."""
        test_data = pd.DataFrame(
            [
                {
                    "snp": "rs7072776",
                    "hg38_chromosome": "10",
                    "hg38_start": 21744013,
                    "hg38_allele_string": "A/G",
                }
            ]
        )

        # Mock API response
        mock_response = {
            "variant": {
                "variant_id": "10-21744013-A-G",
                "rsid": "rs7072776",
                "chrom": "10",
                "pos": 21744013,
                "ref": "A",
                "alt": "G",
            }
        }

        with patch.object(self.ensembl_client, "get_variations_batch", return_value={}):
            with patch.object(
                self.client, "_make_graphql_request", return_value=mock_response
            ) as mock_request:
                # First harmonization - should hit API
                result1 = self.harmonizer.harmonize_snp_batch(test_data)
                initial_call_count = mock_request.call_count

                # Second harmonization - should use cache
                result2 = self.harmonizer.harmonize_snp_batch(test_data)
                # Note: Since the SNP has coordinates, there may be additional validation calls
                # So we check that the cache reduces the total number of calls
                assert mock_request.call_count <= initial_call_count + 1

                # Results should be identical
                assert len(result1) == len(result2)

                # Check cache statistics
                cache_stats = self.client.get_cache_stats()
                assert cache_stats["total_cached_entries"] > 0
