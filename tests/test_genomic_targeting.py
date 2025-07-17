"""
Tests for genomic targeting functionality.

This module tests the genomic targeting flag feature including:
- Configuration validation
- File processing (Excel, CSV, TSV)
- Flag application to gene DataFrames
- Pipeline integration
- Error handling and edge cases
"""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

from custom_panel.core.config_manager import ConfigManager
from custom_panel.sources.genomic_targeting import (
    apply_genomic_targeting_flags,
    convert_to_boolean,
    fetch_genomic_targeting_flags,
    get_genomic_targeting_summary,
    process_targeting_file,
    validate_genomic_targeting_config,
)


class TestGenomicTargetingConfig:
    """Test genomic targeting configuration handling."""

    def test_config_manager_genomic_targeting_methods(self):
        """Test ConfigManager genomic targeting methods."""
        config = {
            "genomic_targeting": {
                "enabled": True,
                "file_path": "test_file.xlsx",
                "gene_column": "gene_symbol",
                "targeting_column": "targeting",
                "default_value": False,
                "allow_missing_file": True,
                "validate_gene_symbols": False,
            }
        }

        config_manager = ConfigManager(config)

        assert config_manager.is_genomic_targeting_enabled() is True
        assert config_manager.get_genomic_targeting_file_path() == "test_file.xlsx"
        assert config_manager.get_genomic_targeting_gene_column() == "gene_symbol"
        assert config_manager.get_genomic_targeting_column() == "targeting"
        assert config_manager.get_genomic_targeting_default_value() is False
        assert config_manager.is_genomic_targeting_missing_file_allowed() is True
        assert config_manager.is_genomic_targeting_validation_enabled() is False

    def test_config_manager_defaults(self):
        """Test ConfigManager defaults when genomic targeting not configured."""
        config = {}
        config_manager = ConfigManager(config)

        assert config_manager.is_genomic_targeting_enabled() is False
        assert (
            config_manager.get_genomic_targeting_file_path()
            == "data/manual/genomic_targeting_flags.xlsx"
        )
        assert config_manager.get_genomic_targeting_gene_column() == "gene_symbol"
        assert config_manager.get_genomic_targeting_column() == "targeting"
        assert config_manager.get_genomic_targeting_default_value() is False
        assert config_manager.is_genomic_targeting_missing_file_allowed() is True
        assert config_manager.is_genomic_targeting_validation_enabled() is False


class TestBooleanConversion:
    """Test boolean conversion functionality."""

    def test_convert_to_boolean_string_values(self):
        """Test conversion of string values to boolean."""
        # True values
        assert convert_to_boolean("true") is True
        assert convert_to_boolean("TRUE") is True
        assert convert_to_boolean("yes") is True
        assert convert_to_boolean("YES") is True
        assert convert_to_boolean("y") is True
        assert convert_to_boolean("Y") is True
        assert convert_to_boolean("1") is True
        assert convert_to_boolean("t") is True
        assert convert_to_boolean("T") is True

        # False values
        assert convert_to_boolean("false") is False
        assert convert_to_boolean("FALSE") is False
        assert convert_to_boolean("no") is False
        assert convert_to_boolean("NO") is False
        assert convert_to_boolean("n") is False
        assert convert_to_boolean("N") is False
        assert convert_to_boolean("0") is False
        assert convert_to_boolean("f") is False
        assert convert_to_boolean("F") is False

    def test_convert_to_boolean_numeric_values(self):
        """Test conversion of numeric values to boolean."""
        assert convert_to_boolean(1) is True
        assert convert_to_boolean(1.0) is True
        assert convert_to_boolean(0) is False
        assert convert_to_boolean(0.0) is False
        assert convert_to_boolean(5) is True  # Non-zero numbers are True

    def test_convert_to_boolean_boolean_values(self):
        """Test conversion of boolean values (should pass through)."""
        assert convert_to_boolean(True) is True
        assert convert_to_boolean(False) is False

    def test_convert_to_boolean_edge_cases(self):
        """Test conversion of edge cases."""
        assert convert_to_boolean(None, default=True) is True
        assert convert_to_boolean(pd.NA, default=False) is False
        assert convert_to_boolean("unknown", default=True) is True
        assert convert_to_boolean("", default=False) is False
        assert convert_to_boolean(" ", default=True) is True  # Whitespace


class TestFileProcessing:
    """Test file processing functionality."""

    def create_test_file(self, file_format: str, data: dict) -> Path:
        """Create a temporary test file with given data."""
        temp_file = tempfile.NamedTemporaryFile(suffix=f".{file_format}", delete=False)
        temp_path = Path(temp_file.name)
        temp_file.close()

        df = pd.DataFrame(data)

        if file_format == "xlsx":
            df.to_excel(temp_path, index=False)
        elif file_format == "csv":
            df.to_csv(temp_path, index=False)
        elif file_format == "tsv":
            df.to_csv(temp_path, sep="\t", index=False)

        return temp_path

    def test_process_targeting_file_excel(self):
        """Test processing Excel targeting file."""
        data = {
            "gene_symbol": ["BRCA1", "BRCA2", "TP53", "APC"],
            "targeting": [True, "yes", 1, "no"],
            "notes": ["Note 1", "Note 2", "Note 3", "Note 4"],
        }

        temp_file = self.create_test_file("xlsx", data)

        try:
            config = {
                "gene_column": "gene_symbol",
                "targeting_column": "targeting",
                "default_value": False,
            }

            result = process_targeting_file(temp_file, config)

            assert len(result) == 4
            assert result["BRCA1"] is True
            assert result["BRCA2"] is True
            assert result["TP53"] is True
            assert result["APC"] is False

        finally:
            temp_file.unlink()

    def test_process_targeting_file_csv(self):
        """Test processing CSV targeting file."""
        data = {
            "gene_symbol": ["MLH1", "MSH2", "CHEK2"],
            "targeting": ["yes", "no", "1"],
        }

        temp_file = self.create_test_file("csv", data)

        try:
            config = {
                "gene_column": "gene_symbol",
                "targeting_column": "targeting",
                "default_value": False,
            }

            result = process_targeting_file(temp_file, config)

            assert len(result) == 3
            assert result["MLH1"] is True
            assert result["MSH2"] is False
            assert result["CHEK2"] is True

        finally:
            temp_file.unlink()

    def test_process_targeting_file_missing_columns(self):
        """Test handling of missing columns in targeting file."""
        data = {"wrong_gene_column": ["BRCA1", "BRCA2"], "targeting": [True, False]}

        temp_file = self.create_test_file("csv", data)

        try:
            config = {
                "gene_column": "gene_symbol",  # This column doesn't exist
                "targeting_column": "targeting",
                "default_value": False,
            }

            with pytest.raises(ValueError, match="Gene column 'gene_symbol' not found"):
                process_targeting_file(temp_file, config)

        finally:
            temp_file.unlink()

    def test_process_targeting_file_empty(self):
        """Test handling of empty targeting file."""
        data = {"gene_symbol": [], "targeting": []}

        temp_file = self.create_test_file("csv", data)

        try:
            config = {
                "gene_column": "gene_symbol",
                "targeting_column": "targeting",
                "default_value": False,
            }

            result = process_targeting_file(temp_file, config)
            assert len(result) == 0

        finally:
            temp_file.unlink()


class TestFlagApplication:
    """Test application of targeting flags to gene DataFrames."""

    def test_apply_genomic_targeting_flags_enabled(self):
        """Test applying targeting flags when feature is enabled."""
        gene_df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "BRCA2", "TP53", "APC"],
                "gene_id": ["ENSG001", "ENSG002", "ENSG003", "ENSG004"],
            }
        )

        targeting_data = {"gene_symbol": ["BRCA1", "TP53"], "targeting": [True, True]}

        temp_file = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
        temp_path = Path(temp_file.name)
        temp_file.close()

        try:
            pd.DataFrame(targeting_data).to_csv(temp_path, index=False)

            config = {
                "genomic_targeting": {
                    "enabled": True,
                    "file_path": str(temp_path),
                    "gene_column": "gene_symbol",
                    "targeting_column": "targeting",
                    "default_value": False,
                    "allow_missing_file": True,
                }
            }

            result = apply_genomic_targeting_flags(gene_df, config)

            assert "genomic_targeting" in result.columns
            assert len(result) == 4

            # Check specific targeting flags
            brca1_row = result[result["approved_symbol"] == "BRCA1"]
            assert brca1_row["genomic_targeting"].iloc[0] is True

            tp53_row = result[result["approved_symbol"] == "TP53"]
            assert tp53_row["genomic_targeting"].iloc[0] is True

            brca2_row = result[result["approved_symbol"] == "BRCA2"]
            assert brca2_row["genomic_targeting"].iloc[0] is False  # Default value

            apc_row = result[result["approved_symbol"] == "APC"]
            assert apc_row["genomic_targeting"].iloc[0] is False  # Default value

        finally:
            temp_path.unlink()

    def test_apply_genomic_targeting_flags_disabled(self):
        """Test behavior when genomic targeting is disabled."""
        gene_df = pd.DataFrame(
            {"approved_symbol": ["BRCA1", "BRCA2"], "gene_id": ["ENSG001", "ENSG002"]}
        )

        config = {"genomic_targeting": {"enabled": False, "default_value": False}}

        result = apply_genomic_targeting_flags(gene_df, config)

        # Should return DataFrame with genomic_targeting column added with default values
        assert "genomic_targeting" in result.columns
        assert all(val is False for val in result["genomic_targeting"])
        assert len(result) == len(gene_df)

    def test_apply_genomic_targeting_flags_missing_approved_symbol(self):
        """Test handling when approved_symbol column is missing."""
        gene_df = pd.DataFrame(
            {
                "gene_symbol": ["BRCA1", "BRCA2"],  # Wrong column name
                "gene_id": ["ENSG001", "ENSG002"],
            }
        )

        config = {"genomic_targeting": {"enabled": True, "default_value": False}}

        result = apply_genomic_targeting_flags(gene_df, config)

        # Should return original DataFrame when approved_symbol column missing
        assert result.equals(gene_df)

    def test_apply_genomic_targeting_flags_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        gene_df = pd.DataFrame()

        config = {"genomic_targeting": {"enabled": True, "default_value": False}}

        result = apply_genomic_targeting_flags(gene_df, config)
        assert result.empty


class TestConfigValidation:
    """Test configuration validation functionality."""

    def test_validate_genomic_targeting_config_valid(self):
        """Test validation of valid configuration."""
        config = {
            "genomic_targeting": {
                "enabled": True,
                "file_path": "test.xlsx",
                "gene_column": "gene_symbol",
                "targeting_column": "targeting",
                "default_value": False,
                "allow_missing_file": True,
                "validate_gene_symbols": False,
            }
        }

        errors = validate_genomic_targeting_config(config)
        assert len(errors) == 0

    def test_validate_genomic_targeting_config_disabled(self):
        """Test validation when feature is disabled."""
        config = {"genomic_targeting": {"enabled": False}}

        errors = validate_genomic_targeting_config(config)
        assert len(errors) == 0  # No validation errors when disabled

    def test_validate_genomic_targeting_config_missing_file_path(self):
        """Test validation with missing file_path."""
        config = {
            "genomic_targeting": {
                "enabled": True
                # Missing file_path
            }
        }

        errors = validate_genomic_targeting_config(config)
        assert len(errors) == 1
        assert "file_path is required" in errors[0]

    def test_validate_genomic_targeting_config_invalid_types(self):
        """Test validation with invalid data types."""
        config = {
            "genomic_targeting": {
                "enabled": "true",  # Should be boolean
                "file_path": 123,  # Should be string
                "gene_column": [],  # Should be string
                "default_value": "false",  # Should be boolean
            }
        }

        errors = validate_genomic_targeting_config(config)
        assert len(errors) >= 4  # Should have multiple type errors


class TestIntegration:
    """Test integration with pipeline and other components."""

    def test_fetch_genomic_targeting_flags_disabled(self):
        """Test fetching flags when feature is disabled."""
        config = {"genomic_targeting": {"enabled": False}}

        result = fetch_genomic_targeting_flags(config)
        assert len(result) == 0

    def test_fetch_genomic_targeting_flags_missing_file(self):
        """Test fetching flags when file is missing but allowed."""
        config = {
            "genomic_targeting": {
                "enabled": True,
                "file_path": "/nonexistent/file.xlsx",
                "allow_missing_file": True,
            }
        }

        result = fetch_genomic_targeting_flags(config)
        assert (
            len(result) == 0
        )  # Should return empty dict when file missing but allowed

    def test_fetch_genomic_targeting_flags_missing_file_not_allowed(self):
        """Test fetching flags when file is missing and not allowed."""
        config = {
            "genomic_targeting": {
                "enabled": True,
                "file_path": "/nonexistent/file.xlsx",
                "allow_missing_file": False,
            }
        }

        with pytest.raises(FileNotFoundError):
            fetch_genomic_targeting_flags(config)

    def test_get_genomic_targeting_summary(self):
        """Test getting configuration summary."""
        config = {
            "genomic_targeting": {
                "enabled": True,
                "file_path": "/nonexistent/file.xlsx",
                "allow_missing_file": True,
            }
        }

        summary = get_genomic_targeting_summary(config)

        assert "enabled" in summary
        assert "file_path" in summary
        assert "validation_errors" in summary
        assert "total_genes" in summary
        assert summary["enabled"] is True


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_process_targeting_file_with_nan_values(self):
        """Test processing file with NaN values."""
        data = {
            "gene_symbol": ["BRCA1", None, "TP53", ""],
            "targeting": [True, "yes", None, False],
        }

        temp_file = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
        temp_path = Path(temp_file.name)
        temp_file.close()

        try:
            pd.DataFrame(data).to_csv(temp_path, index=False)

            config = {
                "gene_column": "gene_symbol",
                "targeting_column": "targeting",
                "default_value": False,
            }

            result = process_targeting_file(temp_path, config)

            # Should only include valid gene symbols with non-null targeting values
            assert "BRCA1" in result
            assert result["BRCA1"] is True
            # TP53 should be excluded because targeting is None
            # Empty string gene should be excluded
            assert (
                len(result) == 1
            )  # Only BRCA1 has both valid gene symbol and targeting value

        finally:
            temp_path.unlink()

    def test_apply_flags_with_duplicate_genes(self):
        """Test applying flags when gene DataFrame has duplicates."""
        gene_df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "BRCA1", "BRCA2"],  # Duplicate BRCA1
                "gene_id": ["ENSG001", "ENSG001", "ENSG002"],
            }
        )

        targeting_data = {"gene_symbol": ["BRCA1"], "targeting": [True]}

        temp_file = tempfile.NamedTemporaryFile(suffix=".csv", delete=False)
        temp_path = Path(temp_file.name)
        temp_file.close()

        try:
            pd.DataFrame(targeting_data).to_csv(temp_path, index=False)

            config = {
                "genomic_targeting": {
                    "enabled": True,
                    "file_path": str(temp_path),
                    "gene_column": "gene_symbol",
                    "targeting_column": "targeting",
                    "default_value": False,
                    "allow_missing_file": True,
                }
            }

            result = apply_genomic_targeting_flags(gene_df, config)

            # Both BRCA1 rows should have targeting=True
            brca1_rows = result[result["approved_symbol"] == "BRCA1"]
            assert all(val is True for val in brca1_rows["genomic_targeting"])

            brca2_rows = result[result["approved_symbol"] == "BRCA2"]
            assert all(val is False for val in brca2_rows["genomic_targeting"])

        finally:
            temp_path.unlink()

    def test_error_handling_behavior(self):
        """Test error handling behavior without mocking."""
        config = {
            "genomic_targeting": {
                "enabled": True,
                "file_path": "/invalid/path/file.xlsx",
                "allow_missing_file": False,
            }
        }

        # Should raise FileNotFoundError when file missing and not allowed
        with pytest.raises(FileNotFoundError):
            fetch_genomic_targeting_flags(config)

        # Test that allowing missing file returns empty dict
        config["genomic_targeting"]["allow_missing_file"] = True
        result = fetch_genomic_targeting_flags(config)
        assert len(result) == 0


class TestHTMLReportIntegration:
    """Test integration with HTML report generation."""

    def test_genomic_targeting_in_html_report_columns(self):
        """Test that genomic_targeting column is included in HTML report."""
        from custom_panel.core.report_generator import ReportGenerator

        # Create test DataFrame with genomic_targeting column
        gene_df = pd.DataFrame(
            {
                "approved_symbol": ["BRCA1", "BRCA2", "TP53"],
                "hgnc_id": ["HGNC:1100", "HGNC:1101", "HGNC:11998"],
                "gene_size": [82000, 84000, 19000],
                "score": [4.2, 4.1, 3.8],
                "include": [True, True, True],
                "genomic_targeting": [True, True, False],
                "source_count": [5, 4, 6],
            }
        )

        report_gen = ReportGenerator()
        table_data, available_columns, default_visible = report_gen._prepare_table_data(
            gene_df
        )

        # Verify genomic_targeting is included
        assert "genomic_targeting" in available_columns
        assert "genomic_targeting" in default_visible

        # Verify data is properly serialized
        assert len(table_data) == 3
        assert table_data[0]["genomic_targeting"] == 1.0  # True -> 1.0
        assert table_data[1]["genomic_targeting"] == 1.0  # True -> 1.0
        assert table_data[2]["genomic_targeting"] == 0.0  # False -> 0.0

