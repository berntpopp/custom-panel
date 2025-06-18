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
        self.source_group_weights = self.scoring_config.get("source_group_weights", {})
        self.thresholds = self.scoring_config.get("thresholds", {})

        # Quality control settings
        self.qc_config = config.get("quality_control", {})

    def create_master_list(
        self, dataframes: list[pd.DataFrame], output_manager: Any = None
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

        # Log source group statistics
        source_stats = {}
        for df in dataframes:
            if not df.empty and "source_group" in df.columns:
                source_group = df["source_group"].iloc[0]
                source_stats[source_group] = len(df)
        logger.info(f"Source group statistics: {source_stats}")

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
            # Choose columns for duplicate detection based on what's available
            if "source_group" in filtered_df.columns:
                # Pre-aggregated data - use source_group
                subset_cols = ["approved_symbol", "source_group"]
            else:
                # Raw data - use source_name
                subset_cols = ["approved_symbol", "source_name"]

            filtered_df = filtered_df.drop_duplicates(subset=subset_cols, keep="first")

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
            f"Scoring configuration - Source group weights: {self.source_group_weights}"
        )
        logger.debug(
            f"Scoring configuration - Max evidence score: {self.scoring_config.get('max_evidence_score', 'No limit')}"
        )

        # Group by approved gene symbol
        grouped = df.groupby("approved_symbol")

        # Initialize lists for aggregated data
        genes = []
        hgnc_ids = []
        final_scores = []
        source_counts = []
        source_names_list = []
        source_details_list = []

        for gene_symbol, group in grouped:
            # Calculate final score using pre-aggregated data
            score = self._calculate_gene_score(group)

            # Aggregate source information (handle both raw and pre-aggregated data)
            if "source_group" in group.columns:
                source_groups = group["source_group"].unique().tolist()
            else:
                source_groups = group["source_name"].unique().tolist()
            source_count = len(source_groups)
            source_names = ";".join(source_groups)

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
            final_scores.append(score)
            source_counts.append(source_count)
            source_names_list.append(source_names)
            source_details_list.append(source_details)

        # Create aggregated DataFrame
        scored_df = pd.DataFrame(
            {
                "approved_symbol": genes,
                "hgnc_id": hgnc_ids,
                "score": final_scores,
                "source_count": source_counts,
                "source_names": source_names_list,
                "source_details": source_details_list,
            }
        )

        # Sort by score (highest first)
        scored_df = scored_df.sort_values("score", ascending=False)

        logger.info(f"Calculated scores for {len(scored_df)} unique genes")

        # Log score distribution summary
        if "score" in scored_df.columns:
            score_stats = scored_df["score"].describe()
            logger.info(
                f"Score distribution - mean: {score_stats['mean']:.2f}, "
                f"median: {score_stats['50%']:.2f}, max: {score_stats['max']:.2f}"
            )

        if "source_count" in scored_df.columns:
            source_dist = scored_df["source_count"].value_counts().sort_index()
            logger.info(f"Source count distribution: {dict(source_dist)}")

        return scored_df

    def _calculate_gene_score(self, group: pd.DataFrame) -> float:
        """
        Calculate final gene score from pre-aggregated source groups.

        Args:
            group: DataFrame group for one gene (pre-aggregated by source group)

        Returns:
            Final weighted score for the gene
        """
        total_score = 0.0

        for _, row in group.iterrows():
            # Handle both raw data (source_name) and pre-aggregated data (source_group)
            if "source_group" in row:
                source_identifier = row["source_group"]
            else:
                source_identifier = row["source_name"]

            base_score = row.get("source_evidence_score", 1.0)
            internal_confidence = row.get("internal_confidence_score", 1.0)

            # Get weight for this source
            weight = self.source_group_weights.get(source_identifier, 1.0)

            # Calculate contribution: base_score * internal_confidence * weight
            score_contribution = base_score * internal_confidence * weight
            total_score += score_contribution

            logger.debug(
                f"Gene score contribution from {source_identifier}: "
                f"{base_score:.2f} * {internal_confidence:.2f} * {weight:.2f} = {score_contribution:.2f}"
            )

        # Apply maximum score limit if configured
        max_score = self.scoring_config.get("max_evidence_score", float("inf"))
        return min(total_score, max_score)

    # _map_source_to_weight_key method is now obsolete and removed

    def _apply_decision_logic(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply decision thresholds to determine gene inclusion."""
        if df.empty:
            return df

        # Get thresholds
        score_threshold = self.thresholds.get("score_threshold", 2.0)
        min_sources = self.thresholds.get("min_sources", 1)

        # Log threshold configuration
        logger.info(
            f"Applying decision thresholds - Score: {score_threshold}, "
            f"Min sources: {min_sources}"
        )

        # Apply decision logic
        df = df.copy()
        df["include"] = (df["score"] >= score_threshold) & (
            df["source_count"] >= min_sources
        )

        # Log decision results
        total_genes = len(df)
        included_genes = df["include"].sum()

        logger.info(f"Decision results: {included_genes}/{total_genes} genes included")

        # Log examples of included/excluded genes for debugging
        if logger.isEnabledFor(logging.DEBUG):
            # Show top included genes
            if included_genes > 0:
                top_genes = df[df["include"]].nlargest(5, "score")
                logger.debug("Top 5 included genes:")
                for _, gene in top_genes.iterrows():
                    logger.debug(
                        f"  {gene['approved_symbol']}: score={gene['score']:.2f}, sources={gene['source_count']}"
                    )

            # Show some excluded genes close to threshold
            close_to_threshold = df[
                ~df["include"]
                & (df["score"] > score_threshold * 0.8)
                & (df["score"] < score_threshold)
            ]
            if len(close_to_threshold) > 0:
                logger.debug("\nGenes close to threshold but excluded:")
                for _, gene in close_to_threshold.head(5).iterrows():
                    logger.debug(
                        f"  {gene['approved_symbol']}: score={gene['score']:.2f}, sources={gene['source_count']}"
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
            "genes_with_inclusion": df["include"].sum()
            if "include" in df.columns
            else 0,
        }

        # Score distributions
        if "score" in df.columns:
            summary["score_stats"] = {
                "mean": df["score"].mean(),
                "median": df["score"].median(),
                "min": df["score"].min(),
                "max": df["score"].max(),
                "std": df["score"].std(),
            }

        # Source count distribution
        if "source_count" in df.columns:
            source_count_dist = df["source_count"].value_counts().sort_index().to_dict()
            summary["source_count_distribution"] = source_count_dist

        # Top scoring genes
        if "score" in df.columns:
            top_genes = df.nlargest(10, "score")[["approved_symbol", "score"]].to_dict(
                "records"
            )
            summary["top_scoring_genes"] = top_genes

        return summary

    def filter_by_category(self, df: pd.DataFrame, category: str) -> pd.DataFrame:
        """
        Filter genes by inclusion.

        Args:
            df: Master DataFrame
            category: Category to filter by (only "include" is supported now)

        Returns:
            Filtered DataFrame
        """
        if df.empty:
            return df

        if category != "include":
            logger.warning(f"Only 'include' category is supported, got: {category}")
            return pd.DataFrame()

        column_name = "include"
        if column_name not in df.columns:
            logger.warning(f"Column {column_name} not found in DataFrame")
            return pd.DataFrame()

        filtered_df = df[df[column_name]].copy()
        logger.info(f"Filtered to {len(filtered_df)} included genes")
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

        # Export included genes
        filtered_df = self.filter_by_category(df, "include")

        if not filtered_df.empty:
            # Export gene list (simple text file)
            gene_list_file = output_dir / "included_genes.txt"
            with open(gene_list_file, "w") as f:
                for gene in filtered_df["approved_symbol"]:
                    f.write(f"{gene}\n")

            exported_files["included_genes"] = str(gene_list_file)
            logger.info(
                f"Exported {len(filtered_df)} included genes to {gene_list_file}"
            )

        return exported_files
