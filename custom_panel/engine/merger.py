"""
Gene panel merger and scoring engine.

This module provides functionality to merge gene panels from multiple sources,
apply scoring logic, and make inclusion decisions.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from ..core.io import merge_panel_dataframes

logger = logging.getLogger(__name__)


class PanelMerger:
    """Gene panel merger and scoring engine."""

    def __init__(self, config: dict[str, Any]):
        """
        Initialize the panel merger.

        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.scoring_config = config.get("scoring", {})
        self.category_weights = self.scoring_config.get("category_weights", {})
        self.thresholds = self.scoring_config.get("thresholds", {})

        # Quality control settings
        self.qc_config = config.get("quality_control", {})

    def create_master_list(
        self, dataframes: list[pd.DataFrame], output_manager=None
    ) -> pd.DataFrame:
        """
        Create master gene list from multiple source DataFrames.

        Args:
            dataframes: List of standardized DataFrames from different sources
            output_manager: Optional output manager for saving intermediate files

        Returns:
            Master DataFrame with merged and scored genes
        """
        if not dataframes:
            logger.warning("No dataframes provided for merging")
            from ..core.io import STANDARD_COLUMNS

            return pd.DataFrame(columns=STANDARD_COLUMNS)

        logger.info(
            f"Starting merger pipeline with {len(dataframes)} source dataframes"
        )

        # Log source statistics
        source_stats = {}
        for df in dataframes:
            if not df.empty and "source_name" in df.columns:
                source_name = df["source_name"].iloc[0]
                source_stats[source_name] = len(df)
        logger.info(f"Source statistics: {source_stats}")

        # Step 1: Merge all dataframes
        merged_df = self._merge_dataframes(dataframes)

        # Save merged data
        if output_manager:
            output_manager.save_merged_data(merged_df, source_stats)

        # Step 2: Apply quality control filters
        filtered_df = self._apply_quality_control(merged_df)

        # Step 3: Aggregate by gene and calculate scores
        scored_df = self._calculate_scores(filtered_df)

        # Save scored data
        if output_manager:
            scoring_summary = self.get_scoring_summary(scored_df)
            output_manager.save_scored_data(scored_df, scoring_summary)

        # Step 4: Apply decision thresholds
        final_df = self._apply_decision_logic(scored_df)

        logger.info(f"Created master list with {len(final_df)} genes")
        return final_df

    def _merge_dataframes(self, dataframes: list[pd.DataFrame]) -> pd.DataFrame:
        """Merge all source dataframes."""
        # Filter out empty dataframes
        valid_dfs = [df for df in dataframes if not df.empty]

        if not valid_dfs:
            return pd.DataFrame()

        # Merge using the IO utility
        merged_df = merge_panel_dataframes(valid_dfs)

        logger.info(f"Merged {len(valid_dfs)} dataframes into {len(merged_df)} records")

        # Log unique genes count
        if "approved_symbol" in merged_df.columns:
            unique_genes = merged_df["approved_symbol"].nunique()
            logger.info(f"Total unique genes before aggregation: {unique_genes}")

        return merged_df

    def _apply_quality_control(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply quality control filters."""
        if df.empty:
            return df

        original_count = len(df)
        filtered_df = df.copy()

        # Remove rows with missing gene symbols
        if "approved_symbol" in filtered_df.columns:
            filtered_df = filtered_df.dropna(subset=["approved_symbol"])
            filtered_df = filtered_df[filtered_df["approved_symbol"].str.strip() != ""]

        # Remove duplicates if configured
        if self.qc_config.get("remove_duplicates", True):
            # Remove exact duplicates based on gene symbol and source
            filtered_df = filtered_df.drop_duplicates(
                subset=["approved_symbol", "source_name"], keep="first"
            )

        removed_count = original_count - len(filtered_df)
        if removed_count > 0:
            logger.info(f"Quality control removed {removed_count} records")

        return filtered_df

    def _calculate_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """Calculate aggregated scores for each gene."""
        if df.empty:
            return df

        logger.info("Calculating gene scores")
        logger.debug(
            f"Scoring configuration - Category weights: {self.category_weights}"
        )
        logger.debug(
            f"Scoring configuration - Max evidence score: {self.scoring_config.get('max_evidence_score', 'No limit')}"
        )

        # Group by approved gene symbol
        grouped = df.groupby("approved_symbol")

        # Initialize lists for aggregated data
        genes = []
        hgnc_ids = []
        germline_scores = []
        somatic_scores = []
        total_scores = []
        source_counts = []
        source_names_list = []
        source_details_list = []

        for gene_symbol, group in grouped:
            # Explicitly filter evidence by category before scoring
            if "category" in group.columns:
                germline_evidence_df = group[group["category"] == "germline"]
                somatic_evidence_df = group[group["category"] == "somatic"]
            else:
                # Backward compatibility: score all evidence for both categories
                # (maintains previous behavior for data without category columns)
                germline_evidence_df = group.copy()
                somatic_evidence_df = group.copy()

            # Calculate scores for filtered evidence data
            germline_score = 0.0
            if not germline_evidence_df.empty:
                germline_score = self._calculate_category_score(
                    germline_evidence_df, "germline"
                )

            somatic_score = 0.0
            if not somatic_evidence_df.empty:
                somatic_score = self._calculate_category_score(
                    somatic_evidence_df, "somatic"
                )

            total_score = germline_score + somatic_score

            # Aggregate source information
            sources = group["source_name"].unique().tolist()
            source_count = len(sources)
            source_names = ";".join(sources)

            # Aggregate source details
            details = group["source_details"].dropna().tolist()
            source_details = ";".join([d for d in details if d])

            # Get HGNC ID (should be the same for all rows of the same gene)
            hgnc_id = (
                group["hgnc_id"].dropna().iloc[0]
                if not group["hgnc_id"].dropna().empty
                else ""
            )

            genes.append(gene_symbol)
            hgnc_ids.append(hgnc_id)
            germline_scores.append(germline_score)
            somatic_scores.append(somatic_score)
            total_scores.append(total_score)
            source_counts.append(source_count)
            source_names_list.append(source_names)
            source_details_list.append(source_details)

        # Create aggregated DataFrame
        scored_df = pd.DataFrame(
            {
                "approved_symbol": genes,
                "hgnc_id": hgnc_ids,
                "germline_score": germline_scores,
                "somatic_score": somatic_scores,
                "total_score": total_scores,
                "source_count": source_counts,
                "source_names": source_names_list,
                "source_details": source_details_list,
            }
        )

        # Sort by total score (highest first)
        scored_df = scored_df.sort_values("total_score", ascending=False)

        logger.info(f"Calculated scores for {len(scored_df)} unique genes")

        # Log score distribution summary
        if "germline_score" in scored_df.columns:
            germline_stats = scored_df["germline_score"].describe()
            logger.info(
                f"Germline score distribution - mean: {germline_stats['mean']:.2f}, "
                f"median: {germline_stats['50%']:.2f}, max: {germline_stats['max']:.2f}"
            )

        if "somatic_score" in scored_df.columns:
            somatic_stats = scored_df["somatic_score"].describe()
            logger.info(
                f"Somatic score distribution - mean: {somatic_stats['mean']:.2f}, "
                f"median: {somatic_stats['50%']:.2f}, max: {somatic_stats['max']:.2f}"
            )

        if "source_count" in scored_df.columns:
            source_dist = scored_df["source_count"].value_counts().sort_index()
            logger.info(f"Source count distribution: {dict(source_dist)}")

        return scored_df

    def _calculate_category_score(self, group: pd.DataFrame, category: str) -> float:
        """
        Calculate weighted score for a specific category (germline/somatic).

        Args:
            group: DataFrame group for one gene
            category: Category name ("germline" or "somatic")

        Returns:
            Weighted score for the category
        """
        category_weights = self.category_weights.get(category, {})
        total_score = 0.0

        for _, row in group.iterrows():
            source_name = row["source_name"]
            evidence_score = row.get("source_evidence_score", 0.0)

            # Map source name to weight key
            weight_key = self._map_source_to_weight_key(source_name)
            weight = category_weights.get(weight_key, 1.0)

            # Calculate weighted score
            weighted_score = evidence_score * weight
            total_score += weighted_score

        # Apply maximum score limit if configured
        max_score = self.scoring_config.get("max_evidence_score", float("inf"))
        return min(total_score, max_score)

    def _map_source_to_weight_key(self, source_name: str) -> str:
        """
        Map a source name to a weight configuration key.

        Args:
            source_name: Source name from the data

        Returns:
            Weight configuration key
        """
        source_name_lower = source_name.lower()

        # Define mapping rules
        if "panelapp" in source_name_lower:
            return "panelapp"
        elif "acmg" in source_name_lower or "incidental" in source_name_lower:
            return "acmg_incidental"
        elif "cosmic" in source_name_lower:
            return "cosmic"
        elif (
            "commercial" in source_name_lower
            or "illumina" in source_name_lower
            or "agilent" in source_name_lower
        ):
            return "commercial_panels"
        elif "hpo" in source_name_lower:
            return "hpo"
        elif "manual" in source_name_lower or "curation" in source_name_lower:
            return "manual_curation"
        else:
            # Default to inhouse_panels for unrecognized sources
            return "inhouse_panels"

    def _apply_decision_logic(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply decision thresholds to determine gene inclusion."""
        if df.empty:
            return df

        # Get thresholds
        germline_threshold = self.thresholds.get("germline_threshold", 2.0)
        somatic_threshold = self.thresholds.get("somatic_threshold", 1.5)
        min_sources = self.thresholds.get("min_sources", 1)

        # Log threshold configuration
        logger.info(
            f"Applying decision thresholds - Germline: {germline_threshold}, "
            f"Somatic: {somatic_threshold}, Min sources: {min_sources}"
        )

        # Apply decision logic
        df = df.copy()
        df["include_germline"] = (df["germline_score"] >= germline_threshold) & (
            df["source_count"] >= min_sources
        )
        df["include_somatic"] = (df["somatic_score"] >= somatic_threshold) & (
            df["source_count"] >= min_sources
        )
        df["include_any"] = df["include_germline"] | df["include_somatic"]

        # Log decision results
        total_genes = len(df)
        germline_genes = df["include_germline"].sum()
        somatic_genes = df["include_somatic"].sum()
        any_genes = df["include_any"].sum()

        logger.info(
            f"Decision results: {germline_genes}/{total_genes} germline, "
            f"{somatic_genes}/{total_genes} somatic, {any_genes}/{total_genes} total"
        )

        # Log examples of included/excluded genes for debugging
        if logger.isEnabledFor(logging.DEBUG):
            # Show top included germline genes
            if germline_genes > 0:
                top_germline = df[df["include_germline"]].nlargest(5, "germline_score")
                logger.debug("Top 5 included germline genes:")
                for _, gene in top_germline.iterrows():
                    logger.debug(
                        f"  {gene['approved_symbol']}: score={gene['germline_score']:.2f}, sources={gene['source_count']}"
                    )

            # Show some excluded genes close to threshold
            close_to_threshold = df[
                ~df["include_germline"]
                & (df["germline_score"] > germline_threshold * 0.8)
                & (df["germline_score"] < germline_threshold)
            ]
            if len(close_to_threshold) > 0:
                logger.debug("\nGenes close to germline threshold but excluded:")
                for _, gene in close_to_threshold.head(5).iterrows():
                    logger.debug(
                        f"  {gene['approved_symbol']}: score={gene['germline_score']:.2f}, sources={gene['source_count']}"
                    )

        return df

    def get_scoring_summary(self, df: pd.DataFrame) -> dict[str, Any]:
        """
        Generate summary statistics for scoring results.

        Args:
            df: Scored DataFrame

        Returns:
            Summary statistics dictionary
        """
        if df.empty:
            return {"total_genes": 0}

        summary = {
            "total_genes": len(df),
            "genes_with_germline_inclusion": df["include_germline"].sum()
            if "include_germline" in df.columns
            else 0,
            "genes_with_somatic_inclusion": df["include_somatic"].sum()
            if "include_somatic" in df.columns
            else 0,
            "genes_with_any_inclusion": df["include_any"].sum()
            if "include_any" in df.columns
            else 0,
        }

        # Score distributions
        if "germline_score" in df.columns:
            summary["germline_score_stats"] = {
                "mean": df["germline_score"].mean(),
                "median": df["germline_score"].median(),
                "min": df["germline_score"].min(),
                "max": df["germline_score"].max(),
                "std": df["germline_score"].std(),
            }

        if "somatic_score" in df.columns:
            summary["somatic_score_stats"] = {
                "mean": df["somatic_score"].mean(),
                "median": df["somatic_score"].median(),
                "min": df["somatic_score"].min(),
                "max": df["somatic_score"].max(),
                "std": df["somatic_score"].std(),
            }

        # Source count distribution
        if "source_count" in df.columns:
            source_count_dist = df["source_count"].value_counts().sort_index().to_dict()
            summary["source_count_distribution"] = source_count_dist

        # Top scoring genes
        if "total_score" in df.columns:
            top_genes = df.nlargest(10, "total_score")[
                ["approved_symbol", "total_score"]
            ].to_dict("records")
            summary["top_scoring_genes"] = top_genes

        return summary

    def filter_by_category(self, df: pd.DataFrame, category: str) -> pd.DataFrame:
        """
        Filter genes by inclusion category.

        Args:
            df: Master DataFrame
            category: Category to filter by ("germline", "somatic", or "any")

        Returns:
            Filtered DataFrame
        """
        if df.empty:
            return df

        column_name = f"include_{category}"
        if column_name not in df.columns:
            logger.warning(f"Column {column_name} not found in DataFrame")
            return pd.DataFrame()

        filtered_df = df[df[column_name]].copy()
        logger.info(f"Filtered to {len(filtered_df)} genes for {category} category")
        return filtered_df

    def export_gene_lists(
        self, df: pd.DataFrame, output_dir: str | Path
    ) -> dict[str, str]:
        """
        Export gene lists for different categories.

        Args:
            df: Master DataFrame
            output_dir: Output directory

        Returns:
            Dictionary mapping category to output file path
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        exported_files = {}

        # Export different categories
        categories = ["germline", "somatic", "any"]

        for category in categories:
            filtered_df = self.filter_by_category(df, category)

            if not filtered_df.empty:
                # Export gene list (simple text file)
                gene_list_file = output_dir / f"{category}_genes.txt"
                with open(gene_list_file, "w") as f:
                    for gene in filtered_df["approved_symbol"]:
                        f.write(f"{gene}\n")

                exported_files[f"{category}_genes"] = str(gene_list_file)
                logger.info(
                    f"Exported {len(filtered_df)} {category} genes to {gene_list_file}"
                )

        return exported_files
