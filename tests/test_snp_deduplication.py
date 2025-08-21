"""Tests for SNP deduplication in reports."""

import pandas as pd


class TestSNPDeduplicationInReports:
    """Test SNP deduplication functionality in report generation."""

    def test_snp_deduplication_in_reports(self):
        """Test that SNP deduplication works properly in report generation."""
        from custom_panel.core.report_generator import ReportGenerator

        # Create test DataFrame with duplicate SNPs (same VCF ID, different categories)
        deduplicated_snp_data = pd.DataFrame({
            "snp": ["12:21178615:T:A", "1:1000000:G:A", "2:2000000:C:T"],
            "rsid": ["rs4149056", "rs123", "rs456"],
            "source": ["Manual_SNPs; PharmGKB", "Manual_SNPs", "PRS_Catalog"],
            "category": ["manual; pharmacogenomics", "manual", "prs"],
            "hg38_chromosome": ["12", "1", "2"],
            "hg38_start": ["21178615", "1000000", "2000000"],
            "hg38_end": ["21178615", "1000000", "2000000"],
        })

        report_gen = ReportGenerator()

        # Test that prepare table data doesn't create duplicates
        table_data = report_gen._prepare_snp_table_data_from_deduplicated(deduplicated_snp_data)

        # Verify all_snps table has no duplicates
        all_snps = table_data["all_snps"]
        assert len(all_snps) == 3, f"Expected 3 unique SNPs, got {len(all_snps)}"

        # Verify no snp_type column in table data
        for snp_entry in all_snps:
            assert "snp_type" not in snp_entry, f"snp_type column found in table data: {snp_entry.keys()}"
            assert "category" in snp_entry, "category column missing from table data"

        # Verify rs4149056 appears only once
        rs4149056_entries = [s for s in all_snps if s["rsid"] == "rs4149056"]
        assert len(rs4149056_entries) == 1, f"rs4149056 should appear once, found {len(rs4149056_entries)} times"

        # Verify merged category is preserved
        rs4149056_entry = rs4149056_entries[0]
        assert rs4149056_entry["category"] == "manual; pharmacogenomics"

    def test_snp_statistics_from_deduplicated(self):
        """Test SNP statistics calculation from deduplicated data."""
        from custom_panel.core.report_generator import ReportGenerator

        # Create test DataFrame with merged categories
        deduplicated_snp_data = pd.DataFrame({
            "snp": ["12:21178615:T:A", "1:1000000:G:A", "2:2000000:C:T"],
            "rsid": ["rs4149056", "rs123", "rs456"],
            "source": ["Manual_SNPs; PharmGKB", "Manual_SNPs", "PRS_Catalog"],
            "category": ["manual; pharmacogenomics", "manual", "prs"],
            "hg38_chromosome": ["12", "1", "2"],
            "hg38_start": ["21178615", "1000000", "2000000"],
            "hg38_end": ["21178615", "1000000", "2000000"],
        })

        report_gen = ReportGenerator()
        stats = report_gen._calculate_snp_statistics_from_deduplicated(deduplicated_snp_data)

        # Verify total count
        assert stats["total_snps"] == 3

        # Verify category counts (should count each category separately)
        category_counts = stats["category_type_counts"]
        assert category_counts["manual"] == 2  # rs4149056 and rs123
        assert category_counts["pharmacogenomics"] == 1  # rs4149056 only
        assert category_counts["prs"] == 1  # rs456 only

        # Verify annotation rate
        assert stats["annotated_snps"] == 3  # All have coordinates
        assert stats["annotation_rate"] == "100.0%"

    def test_no_snp_type_in_deduplicated_data(self):
        """Test that snp_type column is never included in deduplicated data."""
        from custom_panel.engine.pipeline import Pipeline

        # Create test data with snp_type column (simulating old format)
        test_df = pd.DataFrame({
            "snp": ["12:21178615:T:A", "12:21178615:T:A"],  # Same SNP ID
            "rsid": ["rs4149056", "rs4149056"],
            "source": ["Manual_SNPs", "PharmGKB"],
            "category": ["manual", "pharmacogenomics"],
            "snp_type": ["manual_snps", "pharmacogenomics"],  # This should be removed
            "hg38_chromosome": ["12", "12"],
            "hg38_start": ["21178615", "21178615"],
            "hg38_end": ["21178615", "21178615"],
        })

        # Create a mock pipeline to test deduplication
        mock_config = {}
        pipeline = Pipeline(mock_config)

        # Test deduplication removes snp_type column
        deduplicated_df = pipeline._deduplicate_snps_in_pipeline(test_df)

        # Verify snp_type column is not present
        assert "snp_type" not in deduplicated_df.columns, f"snp_type column found in deduplicated data: {deduplicated_df.columns.tolist()}"

        # Verify category column is properly merged
        assert len(deduplicated_df) == 1, "Should have one deduplicated entry"
        assert deduplicated_df.iloc[0]["category"] == "manual; pharmacogenomics"

        # Verify source is properly merged
        assert deduplicated_df.iloc[0]["source"] == "Manual_SNPs; PharmGKB"

    def test_report_table_data_columns(self):
        """Test that report table data only includes expected columns."""
        from custom_panel.core.report_generator import ReportGenerator

        # Create clean deduplicated data
        deduplicated_snp_data = pd.DataFrame({
            "snp": ["12:21178615:T:A"],
            "rsid": ["rs4149056"],
            "source": ["Manual_SNPs; PharmGKB"],
            "category": ["manual; pharmacogenomics"],
            "hg38_chromosome": ["12"],
            "hg38_start": ["21178615"],
            "hg38_end": ["21178615"],
        })

        report_gen = ReportGenerator()
        table_data_list = report_gen._convert_snp_df_to_table_data(deduplicated_snp_data)

        # Verify expected columns are present
        expected_columns = {"snp", "rsid", "source", "category", "hg38_chromosome", "hg38_start", "hg38_end"}

        assert len(table_data_list) == 1
        row = table_data_list[0]

        # Check that all expected columns are present
        for col in expected_columns:
            assert col in row, f"Expected column '{col}' missing from row: {row.keys()}"

        # Check that no unexpected columns are present (especially snp_type)
        for col in row.keys():
            assert col in expected_columns, f"Unexpected column '{col}' found in row data"
