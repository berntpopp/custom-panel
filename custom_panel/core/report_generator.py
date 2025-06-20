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
        self, df: pd.DataFrame, config: dict[str, Any], output_path: Path
    ) -> None:
        """
        Render the HTML report and save it to a file.

        Args:
            df: Final annotated DataFrame
            config: Configuration dictionary
            output_path: Path where the HTML report will be saved

        Raises:
            FileNotFoundError: If the template file is not found
            Exception: If rendering fails
        """
        try:
            template = self.env.get_template("report.html.j2")
            context = self._prepare_context(df, config)
            html_content = template.render(context)

            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(html_content)

            logger.info(f"Generated interactive HTML report: {output_path}")

        except Exception as e:
            logger.error(f"Failed to generate HTML report: {e}")
            raise

    def _prepare_context(
        self, df: pd.DataFrame, config: dict[str, Any]
    ) -> dict[str, Any]:
        """
        Prepare template context from DataFrame and configuration.

        Args:
            df: Final annotated DataFrame
            config: Configuration dictionary

        Returns:
            Dictionary containing template variables
        """
        # Calculate basic statistics
        basic_stats = self._calculate_basic_statistics(df)

        # Calculate source statistics
        source_stats, unique_sources = self._calculate_source_statistics(
            df, basic_stats["total_genes"]
        )
        source_diversity = self._calculate_source_diversity(df, unique_sources)

        # Prepare interactive data
        top_genes = self._get_top_genes(df)
        table_data, available_columns, default_visible = self._prepare_table_data(df)
        chart_data = self._prepare_chart_data(df, source_stats)

        # Combine all context data
        context = {
            "generation_date": datetime.now().strftime("%B %d, %Y at %I:%M %p"),
            **basic_stats,
            **source_diversity,
            "source_stats": source_stats,
            "top_genes": top_genes,
            "table_data": json.dumps(table_data),
            "chart_data": json.dumps(chart_data),
            "available_columns": json.dumps(available_columns),
            "default_visible": json.dumps(default_visible),
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

    def _calculate_source_diversity(
        self, df: pd.DataFrame, unique_sources: set[str]
    ) -> dict[str, Any]:
        """Calculate source diversity metrics."""
        return {
            "total_unique_sources": len(unique_sources),
            "avg_sources_per_gene": f"{safe_column_mean(df, 'source_count'):.1f}",
            "max_sources_per_gene": safe_column_max(df, "source_count", default=0),
        }

    def _calculate_source_statistics(
        self, df: pd.DataFrame, total_genes: int
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
                    "percentage": round((gene_count / total_genes) * 100, 1)
                    if total_genes > 0
                    else 0,
                }

        # Sort sources by gene count and return top 12
        sorted_sources = sorted(
            source_stats.items(), key=lambda x: x[1]["gene_count"], reverse=True
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
                        }
                    )
        return top_genes

    def _prepare_table_data(
        self, df: pd.DataFrame
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
        self, df: pd.DataFrame, columns: list[str]
    ) -> list[dict[str, Any]]:
        """Convert DataFrame to JSON-serializable records."""
        table_data = []

        for _, row in df.iterrows():
            record = self._process_table_row(row, columns)
            table_data.append(record)

        return table_data

    def _process_table_row(
        self, row: pd.Series[Any], columns: list[str]
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
                record["source_names_tooltip"]
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
        self, df: pd.DataFrame, source_stats: list[dict[str, Any]]
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
