"""
Improved HTML report generator for the custom-panel tool.

This module provides functionality to generate interactive HTML reports
from gene panel data using Jinja2 templates with better modularity and DRY principles.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd
else:
    import pandas as pd

from jinja2 import Environment, FileSystemLoader

from .dataframe_utils import (
    extract_numeric_data,
    filter_existing_columns,
    safe_bool_count,
    safe_column_count,
    safe_column_max,
    safe_column_mean,
    safe_float_conversion,
    safe_int_conversion,
    safe_str_conversion,
)

logger = logging.getLogger(__name__)


class ReportGenerator:
    """Generates an interactive HTML report from the final gene panel data."""

    def __init__(self, template_dir: Path | None = None):
        """
        Initialize the report generator.

        Args:
            template_dir: Directory containing Jinja2 templates.
                         If None, uses the default templates directory.
        """
        if template_dir is None:
            template_dir = Path(__file__).parent / "templates"

        self.template_dir = template_dir
        self.env = Environment(
            loader=FileSystemLoader(template_dir),
            autoescape=True,
            trim_blocks=True,
            lstrip_blocks=True,
        )

    def render(
        self,
        df: pd.DataFrame,
        config: dict[str, Any],
        output_path: Path,
        deduplicated_snp_data: pd.DataFrame | None = None,
        regions_data: dict[str, pd.DataFrame] | None = None,
    ) -> None:
        """
        Render the HTML report and save it to a file.

        Args:
            df: Final annotated DataFrame
            config: Configuration dictionary
            output_path: Path where the HTML report will be saved
            deduplicated_snp_data: Optional deduplicated SNP data
            regions_data: Optional regions data dictionary by type

        Raises:
            FileNotFoundError: If the template file is not found
            Exception: If rendering fails
        """
        try:
            template = self.env.get_template("report.html.j2")
            context = self._prepare_context(
                df, config, deduplicated_snp_data, regions_data
            )
            html_content = template.render(context)

            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(html_content)

            logger.info(f"Generated interactive HTML report: {output_path}")

        except Exception as e:
            logger.error(f"Failed to generate HTML report: {e}")
            raise

    def _prepare_context(
        self,
        df: pd.DataFrame,
        config: dict[str, Any],
        deduplicated_snp_data: pd.DataFrame | None = None,
        regions_data: dict[str, pd.DataFrame] | None = None,
    ) -> dict[str, Any]:
        """
        Prepare template context from DataFrame and configuration.

        Args:
            df: Final annotated DataFrame
            config: Configuration dictionary
            deduplicated_snp_data: Optional deduplicated SNP data
            regions_data: Optional regions data dictionary by type

        Returns:
            Dictionary containing template variables
        """
        # Calculate basic statistics
        basic_stats = self._calculate_basic_statistics(df)

        # Calculate gene coverage statistics for panel genes
        gene_coverage_stats = self._calculate_gene_coverage_statistics(df)

        # Calculate source statistics
        source_stats, unique_sources = self._calculate_source_statistics(
            df,
            basic_stats["total_genes"],
        )
        source_diversity = self._calculate_source_diversity(df, unique_sources)

        # Prepare interactive data
        top_genes = self._get_top_genes(df)
        table_data, available_columns, default_visible = self._prepare_table_data(df)
        chart_data = self._prepare_chart_data(df, source_stats)

        # Prepare gene source data for tabbed display
        gene_source_data = self._prepare_gene_source_data(table_data, source_stats)
        gene_source_stats = self._calculate_gene_source_statistics(
            source_stats,
            basic_stats["total_genes"],
        )

        # Add panel genes count for second tab
        panel_genes_count = len([gene for gene in table_data if gene.get("include")])
        gene_source_stats["panel_genes"] = panel_genes_count

        # Prepare SNP data if available
        snp_stats = {}
        snp_table_data = {}
        if deduplicated_snp_data is not None and not deduplicated_snp_data.empty:
            snp_stats = self._calculate_snp_statistics_from_deduplicated(
                deduplicated_snp_data
            )
            snp_table_data = self._prepare_snp_table_data_from_deduplicated(
                deduplicated_snp_data
            )

        # Prepare regions data if available
        regions_stats = {}
        regions_table_data = {}
        if regions_data:
            regions_stats = self._calculate_regions_statistics(regions_data)
            regions_table_data = self._prepare_regions_table_data(regions_data)

        # Combine all context data
        context = {
            "generation_date": datetime.now().strftime("%B %d, %Y at %I:%M %p"),
            **basic_stats,
            **source_diversity,
            "gene_coverage_stats": gene_coverage_stats,
            "source_stats": source_stats,
            "top_genes": top_genes,
            "table_data": json.dumps(table_data),
            "chart_data": json.dumps(chart_data),
            "available_columns": json.dumps(available_columns),
            "default_visible": json.dumps(default_visible),
            # Gene source data
            "gene_source_data": json.dumps(gene_source_data),
            "gene_source_stats": gene_source_stats,
            # SNP data
            "has_snp_data": bool(
                deduplicated_snp_data is not None and not deduplicated_snp_data.empty
            ),
            "snp_stats": snp_stats,
            "snp_table_data": json.dumps(snp_table_data),
            # Regions data
            "has_regions_data": bool(regions_data),
            "regions_stats": regions_stats,
            "regions_table_data": json.dumps(regions_table_data),
        }

        return context

    def _calculate_basic_statistics(self, df: pd.DataFrame) -> dict[str, Any]:
        """Calculate basic summary statistics for the report."""
        return {
            "total_genes": len(df),
            "included_count": safe_bool_count(df, "include"),
            "annotated_count": safe_column_count(df, "chromosome"),
            "mane_select_count": safe_column_count(df, "mane_select_transcript"),
            "mane_clinical_count": safe_column_count(df, "mane_clinical_transcript"),
        }

    def _calculate_gene_coverage_statistics(self, df: pd.DataFrame) -> dict[str, Any]:
        """
        Calculate total coverage for panel genes based on genomic targeting strategy.

        For genes with genomic_targeting = 1: use gene_size
        For genes with genomic_targeting = 0 or missing: use mane_select_coverage
        Only includes genes where include = True
        """
        total_bp = 0
        panel_genes_count = 0
        genomic_targeted_count = 0
        mane_targeted_count = 0

        # Filter to only included genes (panel genes)
        if "include" in df.columns:
            panel_df = df[df["include"]].copy()
        else:
            panel_df = df.copy()

        panel_genes_count = len(panel_df)

        for _, row in panel_df.iterrows():
            # Check genomic targeting flag
            genomic_targeting_value = row.get("genomic_targeting")
            genomic_targeting = (
                genomic_targeting_value == 1
                or genomic_targeting_value is True
                or (
                    isinstance(genomic_targeting_value, str)
                    and genomic_targeting_value.lower() in ["true", "1", "yes"]
                )
            )

            if genomic_targeting:
                # Use gene size for genomically targeted genes
                gene_size = safe_int_conversion(row.get("gene_size"))
                if gene_size is not None:
                    total_bp += gene_size
                    genomic_targeted_count += 1
            else:
                # Use MANE select coverage for transcript-targeted genes
                mane_coverage = safe_int_conversion(row.get("mane_select_coverage"))
                if mane_coverage is not None:
                    total_bp += mane_coverage
                    mane_targeted_count += 1

        # Format coverage in human-readable format
        if total_bp >= 1_000_000:
            total_mb = f"{total_bp / 1_000_000:.1f} MB"
        elif total_bp >= 1_000:
            total_mb = f"{total_bp / 1_000:.1f} KB"
        else:
            total_mb = f"{total_bp} bp"

        return {
            "total_bp": total_bp,
            "total_mb": total_mb,
            "panel_genes_count": panel_genes_count,
            "genomic_targeted_count": genomic_targeted_count,
            "mane_targeted_count": mane_targeted_count,
        }

    def _calculate_source_diversity(
        self,
        df: pd.DataFrame,
        unique_sources: set[str],
    ) -> dict[str, Any]:
        """Calculate source diversity metrics."""
        return {
            "total_unique_sources": len(unique_sources),
            "avg_sources_per_gene": f"{safe_column_mean(df, 'source_count'):.1f}",
            "max_sources_per_gene": safe_column_max(df, "source_count", default=0),
        }

    def _calculate_source_statistics(
        self,
        df: pd.DataFrame,
        total_genes: int,
    ) -> tuple[list[dict[str, Any]], set[str]]:
        """Calculate statistics for data sources."""
        source_stats = {}
        unique_sources = set()
        source_gene_counts: dict[str, set[str]] = {}

        if "source_names" in df.columns:
            # Parse source information from each gene
            for _, row in df.iterrows():
                source_names = safe_str_conversion(row.get("source_names"))
                if source_names:
                    sources = source_names.split(";")
                    gene_symbol = safe_str_conversion(row["approved_symbol"])

                    for source in sources:
                        source = source.strip()
                        if source:
                            unique_sources.add(source)
                            if source not in source_gene_counts:
                                source_gene_counts[source] = set()
                            source_gene_counts[source].add(gene_symbol)

            # Calculate statistics per source
            for source in sorted(unique_sources):
                gene_count = len(source_gene_counts[source])
                source_stats[source] = {
                    "gene_count": gene_count,
                    "percentage": (
                        round((gene_count / total_genes) * 100, 1)
                        if total_genes > 0
                        else 0
                    ),
                }

        # Sort sources by gene count and return top 12
        sorted_sources = sorted(
            source_stats.items(),
            key=lambda x: x[1]["gene_count"],
            reverse=True,
        )

        formatted_sources = [
            {
                "name": source,
                "gene_count": stats["gene_count"],
                "percentage": stats["percentage"],
            }
            for source, stats in sorted_sources[:12]
        ]

        return formatted_sources, unique_sources

    def _get_top_genes(self, df: pd.DataFrame) -> list[dict[str, Any]]:
        """Get top 10 scoring genes."""
        top_genes = []
        if "score" in df.columns and len(df) > 0:
            top_10 = df.nlargest(10, "score")
            for _, row in top_10.iterrows():
                score = safe_float_conversion(row["score"])
                if score is not None:
                    top_genes.append(
                        {
                            "gene": safe_str_conversion(row["approved_symbol"]),
                            "score": score,
                        },
                    )
        return top_genes

    def _prepare_table_data(
        self,
        df: pd.DataFrame,
    ) -> tuple[list[dict[str, Any]], list[str], list[str]]:
        """Prepare data for the interactive DataTable."""
        # Define column sets
        all_potential_columns = [
            "approved_symbol",
            "hgnc_id",
            "gene_size",
            "chromosome",
            "gene_start",
            "gene_end",
            "biotype",
            "gene_description",
            "canonical_transcript_coverage",
            "mane_select_coverage",
            "mane_clinical_coverage",
            "score",
            "include",
            "genomic_targeting",
            "source_count",
            "veto_reasons",
            "inclusion_reason",
        ]

        default_visible_columns = [
            "approved_symbol",
            "gene_size",
            "mane_select_coverage",
            "score",
            "include",
            "genomic_targeting",
            "source_count",
            "inclusion_reason",
        ]

        # Filter to existing columns
        available_columns = filter_existing_columns(df, all_potential_columns)
        default_visible = filter_existing_columns(df, default_visible_columns)

        # Convert DataFrame to records
        table_data = self._convert_df_to_records(df, available_columns)

        return table_data, available_columns, default_visible

    def _convert_df_to_records(
        self,
        df: pd.DataFrame,
        columns: list[str],
    ) -> list[dict[str, Any]]:
        """Convert DataFrame to JSON-serializable records."""
        table_data = []

        for _, row in df.iterrows():
            record = self._process_table_row(row, columns)
            table_data.append(record)

        return table_data

    def _process_table_row(
        self,
        row: pd.Series[Any],
        columns: list[str],
    ) -> dict[str, Any]:
        """Process a single DataFrame row into a table record."""
        record: dict[str, Any] = {}

        # Process each column value
        for col in columns:
            value = row[col]
            record[col] = self._serialize_value(value)

        # Add metadata for tooltips
        record["hgnc_id_tooltip"] = safe_str_conversion(row.get("hgnc_id"))
        record["source_names_tooltip"] = safe_str_conversion(row.get("source_details"))

        # Calculate source count if missing
        if "source_count" not in record or record["source_count"] is None:
            record["source_count"] = self._extract_source_count(
                record["source_names_tooltip"],
            )

        return record

    def _serialize_value(self, value: Any) -> Any:
        """Serialize a value for JSON output."""
        if pd.isna(value):
            return None
        elif isinstance(value, int | float | bool):
            return float(value) if isinstance(value, int | float) else value
        else:
            return str(value)

    def _extract_source_count(self, source_details: str) -> int:
        """Extract source count from source details string."""
        if isinstance(source_details, str) and "sources in" in source_details:
            try:
                count_str = source_details.split(" sources in")[0]
                return int(count_str)
            except Exception:
                pass
        return 1

    def _prepare_chart_data(
        self,
        df: pd.DataFrame,
        source_stats: list[dict[str, Any]],
    ) -> dict[str, Any]:
        """Prepare data for interactive charts."""
        # Extract numeric data efficiently
        numeric_columns = ["score", "gene_size", "source_count"]
        numeric_data = extract_numeric_data(df, numeric_columns)

        # Collect transcript sizes
        transcript_sizes = self._collect_transcript_sizes(df)

        chart_data = {
            "scores": numeric_data["score"],
            "gene_sizes": numeric_data["gene_size"],
            "source_counts": numeric_data["source_count"],
            "transcript_sizes": transcript_sizes,
            "source_labels": [source["name"] for source in source_stats[:10]],
            "source_gene_counts": [
                source["gene_count"] for source in source_stats[:10]
            ],
            "source_percentages": [
                source["percentage"] for source in source_stats[:10]
            ],
        }

        return chart_data

    def _collect_transcript_sizes(self, df: pd.DataFrame) -> list[int]:
        """Collect transcript sizes for distribution chart."""
        transcript_sizes = []

        # Extract from canonical and MANE select coverage columns
        coverage_columns = ["canonical_transcript_coverage", "mane_select_coverage"]

        for col in coverage_columns:
            if col in df.columns:
                for _, row in df.iterrows():
                    coverage = safe_int_conversion(row.get(col))
                    if coverage is not None:
                        transcript_sizes.append(coverage)

        return transcript_sizes

    def _calculate_snp_statistics(
        self,
        snp_data: dict[str, pd.DataFrame],
    ) -> dict[str, Any]:
        """Calculate SNP summary statistics for the report."""
        total_snps = 0
        category_type_counts = {}
        annotated_snps = 0

        for category_type, snp_df in snp_data.items():
            if snp_df.empty:
                continue

            snp_count = len(snp_df)
            total_snps += snp_count
            category_type_counts[category_type] = snp_count

            # Count annotated SNPs (those with genomic coordinates)
            coord_columns = ["hg38_chromosome", "hg38_start", "hg38_end"]
            if all(col in snp_df.columns for col in coord_columns):
                annotated_count = snp_df.dropna(subset=coord_columns).shape[0]
                annotated_snps += annotated_count

        return {
            "total_snps": total_snps,
            "category_types": len(category_type_counts),
            "category_type_counts": category_type_counts,
            "annotated_snps": annotated_snps,
            "annotation_rate": (
                f"{(annotated_snps / total_snps * 100):.1f}%"
                if total_snps > 0
                else "0.0%"
            ),
        }

    def _prepare_snp_table_data(
        self,
        snp_data: dict[str, pd.DataFrame],
    ) -> dict[str, Any]:
        """Prepare SNP data for interactive tables in the HTML report."""
        snp_tables = {}

        # Combine all SNPs for master table
        all_snps = []
        for _, snp_df in snp_data.items():
            if not snp_df.empty:
                snp_df_copy = snp_df.copy()
                # Note: snp_type info already preserved in category column
                all_snps.append(snp_df_copy)

        if all_snps:
            master_snp_df = pd.concat(all_snps, ignore_index=True)
            snp_tables["all_snps"] = self._convert_snp_df_to_table_data(master_snp_df)

            # Individual SNP type tables
            for category_type, snp_df in snp_data.items():
                if not snp_df.empty:
                    snp_tables[f"snps_{category_type}"] = (
                        self._convert_snp_df_to_table_data(snp_df)
                    )

        return snp_tables

    def _convert_snp_df_to_table_data(self, df: pd.DataFrame) -> list[dict[str, Any]]:
        """Convert SNP DataFrame to list of dictionaries for HTML table display."""
        # Select relevant columns for display
        display_columns = ["snp", "rsid", "source", "category"]
        if "hg38_chromosome" in df.columns:
            display_columns.extend(["hg38_chromosome", "hg38_start", "hg38_end"])
        if "hg38_allele_string" in df.columns:
            display_columns.append("hg38_allele_string")
        # Only hg38 coordinates are supported

        # Filter to existing columns
        existing_columns = [col for col in display_columns if col in df.columns]

        # Convert to list of dictionaries
        table_data = []
        for _, row in df.iterrows():
            row_data = {}
            for col in existing_columns:
                value = row.get(col)
                # Convert to string for JSON serialization, handle NaN
                if pd.isna(value):
                    row_data[col] = ""
                else:
                    row_data[col] = str(value)
            table_data.append(row_data)

        return table_data

    def _calculate_snp_statistics_from_deduplicated(
        self,
        deduplicated_snp_data: pd.DataFrame,
    ) -> dict[str, Any]:
        """Calculate SNP summary statistics from deduplicated data."""
        total_snps = len(deduplicated_snp_data)

        # Count SNPs by category (parse from merged categories)
        category_type_counts: dict[str, int] = {}
        if "category" in deduplicated_snp_data.columns:
            for _, row in deduplicated_snp_data.iterrows():
                categories = str(row["category"]).split("; ")
                for cat in categories:
                    cat = cat.strip()
                    if cat:
                        category_type_counts[cat] = category_type_counts.get(cat, 0) + 1

        # Count annotated SNPs (those with genomic coordinates)
        annotated_snps = 0
        coord_columns = ["hg38_chromosome", "hg38_start", "hg38_end"]
        if all(col in deduplicated_snp_data.columns for col in coord_columns):
            annotated_snps = deduplicated_snp_data.dropna(subset=coord_columns).shape[0]

        return {
            "total_snps": total_snps,
            "category_types": len(category_type_counts),
            "category_type_counts": category_type_counts,
            "annotated_snps": annotated_snps,
            "annotation_rate": (
                f"{(annotated_snps / total_snps * 100):.1f}%"
                if total_snps > 0
                else "0.0%"
            ),
        }

    def _prepare_snp_table_data_from_deduplicated(
        self,
        deduplicated_snp_data: pd.DataFrame,
    ) -> dict[str, Any]:
        """Prepare SNP table data from deduplicated data."""
        # Convert the deduplicated DataFrame directly to table data
        all_snps_data = self._convert_snp_df_to_table_data(deduplicated_snp_data)

        # Create separate views for each category
        snp_tables = {"all_snps": all_snps_data}

        # Group by individual categories for separate tabs
        if "category" in deduplicated_snp_data.columns:
            category_groups: dict[str, list[dict[str, Any]]] = {}
            for _, row in deduplicated_snp_data.iterrows():
                categories = str(row["category"]).split("; ")
                for cat in categories:
                    cat = cat.strip()
                    if cat:
                        if cat not in category_groups:
                            category_groups[cat] = []
                        # Convert row to dict for table data
                        row_data = {}
                        for col in ["snp", "rsid", "source", "category"]:
                            if col in row:
                                value = row[col]
                                row_data[col] = "" if pd.isna(value) else str(value)

                        # Add coordinate columns if available
                        for col in ["hg38_chromosome", "hg38_start", "hg38_end"]:
                            if col in row:
                                value = row[col]
                                row_data[col] = "" if pd.isna(value) else str(value)

                        category_groups[cat].append(row_data)

            # Add category-specific tables
            for cat, rows in category_groups.items():
                snp_tables[f"snps_{cat}"] = rows

        return snp_tables

    def _calculate_regions_statistics(
        self,
        regions_data: dict[str, pd.DataFrame],
    ) -> dict[str, Any]:
        """Calculate regions summary statistics for the report."""
        total_regions = 0
        region_type_counts = {}
        total_bp = 0

        for region_type, region_df in regions_data.items():
            if region_df.empty:
                continue

            region_count = len(region_df)
            total_regions += region_count
            region_type_counts[region_type] = region_count

            # Calculate total base pairs covered using coverage column if available
            if "coverage" in region_df.columns:
                for _, row in region_df.iterrows():
                    try:
                        coverage = int(row["coverage"])
                        total_bp += coverage
                    except (ValueError, TypeError):
                        continue
            elif "start" in region_df.columns and "end" in region_df.columns:
                # Fallback to start/end calculation
                for _, row in region_df.iterrows():
                    try:
                        start = int(row["start"])
                        end = int(row["end"])
                        total_bp += end - start + 1
                    except (ValueError, TypeError):
                        continue

        # Format coverage in human-readable format
        if total_bp >= 1_000_000:
            total_mb = f"{total_bp / 1_000_000:.1f} MB"
        elif total_bp >= 1_000:
            total_mb = f"{total_bp / 1_000:.1f} KB"
        else:
            total_mb = f"{total_bp} bp"

        return {
            "total_regions": total_regions,
            "region_types": len(region_type_counts),
            "region_type_counts": region_type_counts,
            "total_bp": total_bp,
            "total_mb": total_mb,
        }

    def _prepare_regions_table_data(
        self,
        regions_data: dict[str, pd.DataFrame],
    ) -> dict[str, Any]:
        """Prepare regions data for interactive tables in the HTML report."""
        regions_tables = {}

        # Combine all regions for master table
        all_regions = []
        for region_type, region_df in regions_data.items():
            if not region_df.empty:
                region_df_copy = region_df.copy()
                region_df_copy["region_type"] = region_type
                all_regions.append(region_df_copy)

        if all_regions:
            master_regions_df = pd.concat(all_regions, ignore_index=True)
            regions_tables["all_regions"] = self._convert_regions_df_to_table_data(
                master_regions_df,
            )

        # Individual region type tables
        for region_type, region_df in regions_data.items():
            if not region_df.empty:
                regions_tables[region_type] = self._convert_regions_df_to_table_data(
                    region_df,
                )

        return regions_tables

    def _convert_regions_df_to_table_data(
        self,
        df: pd.DataFrame,
    ) -> list[dict[str, Any]]:
        """Convert regions DataFrame to list of dictionaries for HTML table display."""
        # Select relevant columns for display in optimal order
        display_columns = [
            "region_name",
            "chromosome",
            "start",
            "end",
            "coverage",
            "region_type",
            "source_type",
            "comment",
        ]

        # Filter to existing columns
        existing_columns = [col for col in display_columns if col in df.columns]

        # Convert to list of dictionaries
        table_data = []
        for _, row in df.iterrows():
            row_data = {}
            for col in existing_columns:
                value = row.get(col)
                # Convert to string for JSON serialization, handle NaN
                if pd.isna(value):
                    row_data[col] = ""
                else:
                    row_data[col] = str(value)
            table_data.append(row_data)

        return table_data

    def _prepare_gene_source_data(
        self,
        table_data: list[dict[str, Any]],
        source_stats: list[dict[str, Any]],
    ) -> dict[str, Any]:
        """
        Prepare gene data organized by source for tabbed display.
        Reuses existing processed table data for efficiency.

        Args:
            table_data: Already processed gene table data
            source_stats: List of source statistics from _calculate_source_statistics

        Returns:
            Dictionary with gene data organized by source
        """
        gene_source_data = {}

        # Add "All Genes" table (reuse existing processed data)
        gene_source_data["all_genes"] = table_data

        # Add "Panel Genes" table (only included genes - second tab)
        panel_genes = []
        for gene_record in table_data:
            if gene_record.get(
                "include",
            ):  # Only genes marked for inclusion (truthy values)
                panel_genes.append(gene_record)

        if panel_genes:
            gene_source_data["panel_genes"] = panel_genes

        # Create source-specific gene lists by filtering existing data
        for source_stat in source_stats:
            source_name = source_stat["name"]
            source_key = (
                f"genes_{source_name.lower().replace(' ', '_').replace('-', '_')}"
            )

            # Filter genes that have this source
            source_genes = []
            for gene_record in table_data:
                gene_sources = gene_record.get("source_names_tooltip", "")
                if gene_sources and source_name in gene_sources:
                    source_genes.append(gene_record)

            # Store source-specific gene data
            if source_genes:
                gene_source_data[source_key] = source_genes

        return gene_source_data

    def _calculate_gene_source_statistics(
        self,
        source_stats: list[dict[str, Any]],
        total_genes: int,
    ) -> dict[str, Any]:
        """
        Calculate statistics for gene sources for tab display.
        Reuses existing source statistics for efficiency.

        Args:
            source_stats: List of source statistics
            total_genes: Total number of genes

        Returns:
            Dictionary with gene source statistics
        """
        gene_source_stats = {}

        # Add "All Genes" count
        gene_source_stats["all_genes"] = total_genes

        # Add individual source counts
        for source_stat in source_stats:
            source_name = source_stat["name"]
            source_key = (
                f"genes_{source_name.lower().replace(' ', '_').replace('-', '_')}"
            )
            gene_source_stats[source_key] = source_stat["gene_count"]

        return gene_source_stats
