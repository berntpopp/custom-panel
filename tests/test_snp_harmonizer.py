"""
Tests for SNP harmonizer functionality.
"""

from unittest.mock import Mock, patch

import pandas as pd
import pytest

from custom_panel.core.ensembl_client import EnsemblClient
from custom_panel.core.snp_harmonizer import SNPHarmonizer


class TestSNPHarmonizer:
    """Test SNP harmonizer functionality."""

    def test_init(self):
        """Test harmonizer initialization."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client, config={"ensembl_batch_size": 50})
        
        assert harmonizer.ensembl_client == client
        assert harmonizer.batch_size == 50
        assert harmonizer.stats["total_processed"] == 0
        assert harmonizer.stats["coordinates_resolved"] == 0

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        empty_df = pd.DataFrame()
        result = harmonizer.harmonize_snp_batch(empty_df)
        
        assert result.empty
        assert harmonizer.stats["total_processed"] == 0

    @patch.object(EnsemblClient, "get_variations_batch")
    def test_successful_ensembl_lookup(self, mock_variations):
        """Test successful rsID resolution via Ensembl."""
        # Mock Ensembl response for normal rsID
        mock_variations.return_value = {
            "rs429358": {
                "mappings": [{
                    "assembly_name": "GRCh38",
                    "seq_region_name": "19",
                    "start": 44908684,
                    "end": 44908684,
                    "strand": 1,
                    "allele_string": "T/C"
                }]
            }
        }
        
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        test_df = pd.DataFrame({
            "snp": ["rs429358"],
            "source": ["test"],
            "category": ["normal"]
        })
        
        result = harmonizer.harmonize_snp_batch(test_df)
        
        assert len(result) == 1
        assert result.iloc[0]["hg38_chromosome"] == "19"
        assert result.iloc[0]["hg38_start"] == 44908684
        assert result.iloc[0]["rsid"] == "rs429358"
        assert harmonizer.stats["coordinates_resolved"] == 1

    @patch.object(EnsemblClient, "get_variations_batch")
    def test_merged_rsid_synonym_mapping(self, mock_variations):
        """Test handling of merged rsIDs via synonym mapping."""
        # Mock Ensembl response where merged rsID is returned as canonical with synonyms
        mock_variations.return_value = {
            "rs4688963": {  # Canonical rsID
                "synonyms": ["rs386594666", "rs58163327"],  # Original merged rsID in synonyms
                "mappings": [{
                    "assembly_name": "GRCh38",
                    "seq_region_name": "4",
                    "start": 5748177,
                    "end": 5748177,
                    "strand": 1,
                    "allele_string": "T/C"
                }]
            }
        }
        
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        test_df = pd.DataFrame({
            "snp": ["rs386594666"],  # Merged rsID
            "source": ["IDT"],
            "category": ["identity"]
        })
        
        result = harmonizer.harmonize_snp_batch(test_df)
        
        assert len(result) == 1
        assert result.iloc[0]["hg38_chromosome"] == "4"
        assert result.iloc[0]["hg38_start"] == 5748177
        assert result.iloc[0]["rsid"] == "rs386594666"  # Original rsID preserved
        assert harmonizer.stats["coordinates_resolved"] == 1

    @patch.object(EnsemblClient, "get_variations_batch")
    def test_source_coordinate_fallback(self, mock_variations):
        """Test fallback to source coordinates when Ensembl fails."""
        # Mock empty Ensembl response (rsID not found)
        mock_variations.return_value = {}
        
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        test_df = pd.DataFrame({
            "snp": ["rs2492591692"],
            "location": ["NC_000010.11:94949144"],  # PharmGKB format
            "source": ["PharmGKB"],
            "category": ["pharmacogenomics"]
        })
        
        result = harmonizer.harmonize_snp_batch(test_df)
        
        assert len(result) == 1
        assert result.iloc[0]["hg38_chromosome"] == "10"
        assert result.iloc[0]["hg38_start"] == 94949144
        assert result.iloc[0]["rsid"] == "rs2492591692"
        assert harmonizer.stats["coordinates_resolved"] == 1

    @patch.object(EnsemblClient, "get_variations_batch")
    def test_combined_mechanisms(self, mock_variations):
        """Test combined merged rsID + source coordinate fallback."""
        # Mock response with one canonical rsID and one missing
        mock_variations.return_value = {
            "rs4688963": {  # Canonical for merged rs386594666
                "synonyms": ["rs386594666"],
                "mappings": [{
                    "assembly_name": "GRCh38",
                    "seq_region_name": "4",
                    "start": 5748177,
                    "end": 5748177,
                    "strand": 1,
                    "allele_string": "T/C"
                }]
            }
            # rs2492591692 missing from Ensembl
        }
        
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        test_df = pd.DataFrame({
            "snp": ["rs386594666", "rs2492591692"],
            "location": [None, "NC_000010.11:94949144"],
            "source": ["IDT", "PharmGKB"],
            "category": ["identity", "pharmacogenomics"]
        })
        
        result = harmonizer.harmonize_snp_batch(test_df)
        
        assert len(result) == 2
        # First SNP: resolved via synonym mapping
        assert result.iloc[0]["hg38_chromosome"] == "4"
        assert result.iloc[0]["hg38_start"] == 5748177
        # Second SNP: resolved via source coordinates
        assert result.iloc[1]["hg38_chromosome"] == "10"
        assert result.iloc[1]["hg38_start"] == 94949144
        assert harmonizer.stats["coordinates_resolved"] == 2


class TestSourceCoordinateExtraction:
    """Test source coordinate extraction functionality."""

    def test_extract_source_coordinates_with_location_column(self):
        """Test extracting coordinates from location column."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        snp_df = pd.DataFrame({
            "snp": ["rs123", "rs456"],
            "location": ["NC_000010.11:94949144", "chr5:12345678"],
            "source": ["PharmGKB", "test"]
        })
        
        result = harmonizer._extract_source_coordinates(["rs123", "rs456"], snp_df)
        
        assert "rs123" in result
        assert result["rs123"]["chromosome"] == "10"
        assert result["rs123"]["start"] == 94949144
        
        assert "rs456" in result
        assert result["rs456"]["chromosome"] == "5"
        assert result["rs456"]["start"] == 12345678

    def test_extract_source_coordinates_missing_rsid(self):
        """Test handling of missing rsID in source data."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        snp_df = pd.DataFrame({
            "snp": ["rs123"],
            "location": ["NC_000010.11:94949144"]
        })
        
        result = harmonizer._extract_source_coordinates(["rs999"], snp_df)
        
        assert result["rs999"] is None

    def test_extract_source_coordinates_no_location(self):
        """Test handling of missing location data."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        snp_df = pd.DataFrame({
            "snp": ["rs123"],
            "source": ["test"]
            # No location column
        })
        
        result = harmonizer._extract_source_coordinates(["rs123"], snp_df)
        
        assert result["rs123"] is None


class TestGenomicLocationParsing:
    """Test genomic location parsing functionality."""

    def test_parse_ncbi_refseq_format(self):
        """Test parsing NCBI RefSeq coordinates."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        # Standard NCBI RefSeq format
        result = harmonizer._parse_genomic_location("NC_000010.11:94949144")
        
        assert result is not None
        assert result["chromosome"] == "10"
        assert result["start"] == 94949144
        assert result["end"] == 94949144
        assert result["assembly"] == "GRCh38"

    def test_parse_chromosome_format(self):
        """Test parsing direct chromosome coordinates."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        # Direct chromosome format
        result = harmonizer._parse_genomic_location("chr10:94949144")
        
        assert result is not None
        assert result["chromosome"] == "10"
        assert result["start"] == 94949144
        assert result["end"] == 94949144

    def test_parse_chromosome_range_format(self):
        """Test parsing chromosome range coordinates."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        # Chromosome range format
        result = harmonizer._parse_genomic_location("chr10:94949144-94949150")
        
        assert result is not None
        assert result["chromosome"] == "10"
        assert result["start"] == 94949144
        assert result["end"] == 94949150

    def test_parse_invalid_format(self):
        """Test handling of invalid coordinate formats."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        # Invalid formats
        assert harmonizer._parse_genomic_location("") is None
        assert harmonizer._parse_genomic_location("invalid") is None
        assert harmonizer._parse_genomic_location("chr10") is None
        assert harmonizer._parse_genomic_location("NC_invalid:pos") is None

    def test_parse_ncbi_refseq_edge_cases(self):
        """Test NCBI RefSeq parsing edge cases."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        # Different chromosome numbers
        result1 = harmonizer._parse_genomic_location("NC_000001.11:12345")
        assert result1["chromosome"] == "1"
        
        result2 = harmonizer._parse_genomic_location("NC_000023.11:67890")  # X chromosome
        assert result2["chromosome"] == "23"

    def test_parse_chromosome_without_chr_prefix(self):
        """Test parsing coordinates without chr prefix."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        result = harmonizer._parse_genomic_location("10:94949144")
        
        assert result is not None
        assert result["chromosome"] == "10"
        assert result["start"] == 94949144


class TestHarmonizerStatistics:
    """Test harmonizer statistics functionality."""

    def test_get_stats(self):
        """Test getting harmonizer statistics."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        stats = harmonizer.get_stats()
        
        assert "total_processed" in stats
        assert "coordinates_resolved" in stats
        assert "errors" in stats
        assert stats["total_processed"] == 0

    def test_reset_stats(self):
        """Test resetting harmonizer statistics."""
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        # Modify stats
        harmonizer.stats["total_processed"] = 10
        harmonizer.stats["coordinates_resolved"] = 8
        
        # Reset
        harmonizer.reset_stats()
        
        assert harmonizer.stats["total_processed"] == 0
        assert harmonizer.stats["coordinates_resolved"] == 0
        assert harmonizer.stats["errors"] == 0


class TestVCFIDGeneration:
    """Test VCF ID generation functionality."""

    @patch.object(EnsemblClient, "get_variations_batch")
    def test_vcf_id_generation(self, mock_variations):
        """Test VCF format ID generation."""
        mock_variations.return_value = {
            "rs429358": {
                "mappings": [{
                    "assembly_name": "GRCh38",
                    "seq_region_name": "19",
                    "start": 44908684,
                    "end": 44908684,
                    "strand": 1,
                    "allele_string": "T/C"
                }]
            }
        }
        
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        test_df = pd.DataFrame({
            "snp": ["rs429358"],
            "source": ["test"]
        })
        
        result = harmonizer.harmonize_snp_batch(test_df)
        
        # Should generate VCF format ID
        assert result.iloc[0]["snp"] == "19:44908684:T:C"
        assert result.iloc[0]["rsid"] == "rs429358"  # Original rsID preserved

    @patch.object(EnsemblClient, "get_variations_batch")
    def test_vcf_id_generation_with_source_coordinates(self, mock_variations):
        """Test VCF ID generation with source coordinates (no alleles)."""
        mock_variations.return_value = {}  # Empty response
        
        client = EnsemblClient()
        harmonizer = SNPHarmonizer(client)
        
        test_df = pd.DataFrame({
            "snp": ["rs2492591692"],
            "location": ["NC_000010.11:94949144"],
            "source": ["PharmGKB"]
        })
        
        result = harmonizer.harmonize_snp_batch(test_df)
        
        # Should generate VCF format ID with N/N alleles
        assert result.iloc[0]["snp"] == "10:94949144:N:N"
        assert result.iloc[0]["rsid"] == "rs2492591692"