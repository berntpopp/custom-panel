"""
Configuration manager for consistent config access throughout the application.

This module provides a centralized way to access configuration values
with type safety and default handling.
"""

from pathlib import Path
from typing import Any, Optional

import yaml


class ConfigManager:
    """Manages configuration access with type safety and defaults."""

    def __init__(self, config: dict[str, Any]):
        """
        Initialize with configuration dictionary.

        Args:
            config: Configuration dictionary
        """
        self.config = config

    def get_nested(self, *keys: str, default: Any = None) -> Any:
        """
        Safely get nested configuration values.

        Args:
            *keys: Sequence of keys to traverse
            default: Default value if key path doesn't exist

        Returns:
            Configuration value or default

        Example:
            config.get_nested("data_sources", "PanelApp", "enabled", default=True)
        """
        result = self.config
        for key in keys:
            if isinstance(result, dict) and key in result:
                result = result[key]
            else:
                return default
        return result

    def get_source_config(self, source_name: str) -> dict[str, Any]:
        """
        Get configuration for a specific data source.

        Args:
            source_name: Name of the data source

        Returns:
            Source configuration dictionary
        """
        return self.get_nested("data_sources", source_name, default={})

    def is_source_enabled(self, source_name: str) -> bool:
        """
        Check if a data source is enabled.

        Args:
            source_name: Name of the data source

        Returns:
            True if source is enabled, False otherwise
        """
        return self.get_nested("data_sources", source_name, "enabled", default=True)

    def get_output_config(self) -> dict[str, Any]:
        """
        Get output configuration section.

        Returns:
            Output configuration dictionary
        """
        return self.get_nested("output", default={})

    def get_output_formats(self) -> list[str]:
        """
        Get list of enabled output formats.

        Returns:
            List of output format strings
        """
        return self.get_nested("output", "formats", default=["excel", "csv", "parquet"])

    def get_output_dir(self) -> str:
        """
        Get output directory path.

        Returns:
            Output directory path string
        """
        return self.get_nested("general", "output_dir", default="results")

    def get_bed_config(self) -> dict[str, Any]:
        """
        Get BED file configuration.

        Returns:
            BED configuration dictionary
        """
        return self.get_nested("output", "bed_files", default={})

    def is_bed_enabled(self, bed_type: str = "germline") -> bool:
        """
        Check if BED file generation is enabled for a specific type.

        Args:
            bed_type: Type of BED file ("germline", "exons", "complete_panel", "genes_all", etc.)

        Returns:
            True if BED generation is enabled
        """
        if bed_type == "exons":
            return self.get_nested(
                "output", "bed_files", "exons", "enabled", default=False
            )
        # Support for new master BED file types
        if bed_type in [
            "complete_panel",
            "complete_panel_exons",
            "complete_panel_genes",
            "genes_all",
            "genes_included",
            "snps_all",
            "regions_all",
        ]:
            return self.get_nested("output", "bed_files", bed_type, default=False)
        # Support for individual category files
        if bed_type == "individual_categories":
            return self.get_nested(
                "output", "bed_files", "individual_categories", default=True
            )
        return self.get_nested("output", "bed_files", bed_type, default=False)

    def get_bed_padding(self) -> int:
        """
        Get BED file padding value.

        Returns:
            Padding value in base pairs
        """
        return self.get_nested("output", "bed_files", "padding", default=0)

    def get_html_config(self) -> dict[str, Any]:
        """
        Get HTML report configuration.

        Returns:
            HTML configuration dictionary
        """
        return self.get_nested("output", "html_report", default={})

    def is_html_enabled(self) -> bool:
        """
        Check if HTML report generation is enabled.

        Returns:
            True if HTML report is enabled
        """
        return self.get_nested("output", "html_report", "enabled", default=True)

    def get_scoring_config(self) -> dict[str, Any]:
        """
        Get scoring configuration section.

        Returns:
            Scoring configuration dictionary
        """
        return self.get_nested("scoring", default={})

    def get_score_threshold(self) -> Optional[float]:
        """
        Get score threshold for gene inclusion.

        Returns:
            Score threshold or None if not set
        """
        return self.get_nested("scoring", "thresholds", "score_threshold")

    def get_min_sources(self) -> Optional[int]:
        """
        Get minimum number of sources required.

        Returns:
            Minimum sources count or None if not set
        """
        return self.get_nested("scoring", "thresholds", "min_sources")

    def get_cache_config(self) -> dict[str, Any]:
        """
        Get cache configuration.

        Returns:
            Cache configuration dictionary
        """
        return self.get_nested("cache", default={})

    def get_cache_ttl(self) -> float:
        """
        Get cache TTL in days.

        Returns:
            Cache TTL in days, defaults to 30.0
        """
        return self.get_nested("cache", "ttl_days", default=30.0)

    def get_cache_dir(self) -> str:
        """
        Get cache directory path.

        Returns:
            Cache directory path string
        """
        return self.get_nested("cache", "cache_dir", default=".cache")

    def get_intermediate_config(self) -> dict[str, Any]:
        """
        Get intermediate files configuration.

        Returns:
            Intermediate files configuration dictionary
        """
        return self.get_nested("output", "intermediate_files", default={})

    def is_intermediate_enabled(self) -> bool:
        """
        Check if intermediate file saving is enabled.

        Returns:
            True if intermediate files should be saved
        """
        return self.get_nested("output", "intermediate_files", "enabled", default=False)

    def get_intermediate_format(self) -> str:
        """
        Get format for intermediate files.

        Returns:
            Intermediate file format string
        """
        return self.get_nested(
            "output", "intermediate_files", "format", default="excel"
        )

    def is_structured_output_enabled(self) -> bool:
        """
        Check if structured output directories are enabled.

        Returns:
            True if structured output should be used
        """
        return self.get_nested(
            "directory_structure", "use_structured_output", default=True
        )

    def get_log_level(self) -> str:
        """
        Get logging level.

        Returns:
            Log level string, defaults to "INFO"
        """
        return self.get_nested("logging", "level", default="INFO")

    def is_file_logging_enabled(self) -> bool:
        """
        Check if file logging is enabled.

        Returns:
            True if logs should be written to files
        """
        return self.get_nested("output", "file_logging", "enabled", default=False)

    def get_source_evidence_score(self, source_name: str) -> float:
        """
        Get evidence score for a specific source.

        Args:
            source_name: Name of the data source

        Returns:
            Evidence score, defaults to 1.0
        """
        return self.get_nested(
            "data_sources", source_name, "evidence_score", default=1.0
        )

    def get_source_category(self, source_name: str) -> str:
        """
        Get category for a specific source.

        Args:
            source_name: Name of the data source

        Returns:
            Source category, defaults to "germline"
        """
        return self.get_nested(
            "data_sources", source_name, "category", default="germline"
        )

    def is_source_group(self, source_name: str) -> bool:
        """
        Check if a source is configured as a source group.

        Args:
            source_name: Name of the data source

        Returns:
            True if source is a source group
        """
        return self.get_nested(
            "data_sources", source_name, "source_group", default=False
        )

    def get_source_normalization(self, source_name: str) -> dict[str, Any]:
        """
        Get normalization configuration for a source.

        Args:
            source_name: Name of the data source

        Returns:
            Normalization configuration dictionary
        """
        return self.get_nested("data_sources", source_name, "normalization", default={})

    def get_exon_padding(self) -> int:
        """
        Get exon padding for BED file generation.

        Returns:
            Exon padding in base pairs, defaults to 10
        """
        return self.get_nested(
            "output", "bed_files", "exons", "exon_padding", default=10
        )

    def get_transcript_types_config(self) -> dict[str, bool]:
        """
        Get configuration for which transcript types to include in exon BED files.

        Returns:
            Dictionary mapping transcript types to enabled status
        """
        exon_config = self.get_nested("output", "bed_files", "exons", default={})
        return {
            "canonical": exon_config.get("canonical_transcript", True),
            "mane_select": exon_config.get("mane_select_transcript", True),
            "mane_clinical": exon_config.get("mane_clinical_transcript", False),
        }

    def override_with_cli_args(self, **kwargs: Any) -> None:
        """
        Override configuration with command line arguments.

        Args:
            **kwargs: Keyword arguments from CLI
        """
        if kwargs.get("output_dir"):
            self._set_nested("general", "output_dir", kwargs["output_dir"])

        if kwargs.get("score_threshold") is not None:
            self._set_nested(
                "scoring", "thresholds", "score_threshold", kwargs["score_threshold"]
            )

        if kwargs.get("save_intermediate"):
            self._set_nested("output", "intermediate_files", "enabled", True)

        if kwargs.get("intermediate_format"):
            self._set_nested(
                "output", "intermediate_files", "format", kwargs["intermediate_format"]
            )

        if kwargs.get("log_to_file"):
            self._set_nested("output", "file_logging", "enabled", True)

        if kwargs.get("structured_output") is False:
            self._set_nested("directory_structure", "use_structured_output", False)

    def _set_nested(self, *keys_and_value: Any) -> None:
        """
        Set a nested configuration value.

        Args:
            *keys_and_value: Keys to traverse and final value to set
        """
        *keys, value = keys_and_value
        target = self.config

        # Navigate to the parent of the target key
        for key in keys[:-1]:
            if key not in target:
                target[key] = {}
            target = target[key]

        # Set the final value
        target[keys[-1]] = value

    @classmethod
    def from_files(
        cls,
        default_path: Path,
        override_path: Optional[Path] = None,
        local_path: Optional[Path] = None,
    ) -> "ConfigManager":
        """Load configuration from files and create a ConfigManager instance."""
        if not default_path.exists():
            raise FileNotFoundError(
                "Default configuration file not found. Installation may be corrupted."
            )

        with open(default_path) as f:
            config = yaml.safe_load(f) or {}

        if override_path and override_path.exists():
            with open(override_path) as f:
                override_config = yaml.safe_load(f) or {}
            config = cls._merge_configs(config, override_config)

        if local_path and local_path.exists():
            with open(local_path) as f:
                local_config = yaml.safe_load(f) or {}
            config = cls._merge_configs(config, local_config)

        return cls(config)

    @staticmethod
    def _merge_configs(
        base_config: dict[str, Any], override_config: dict[str, Any]
    ) -> dict[str, Any]:
        """Recursively merge override configuration into base configuration."""
        import copy

        result = copy.deepcopy(base_config)
        for key, value in override_config.items():
            if (
                key in result
                and isinstance(result.get(key), dict)
                and isinstance(value, dict)
            ):
                result[key] = ConfigManager._merge_configs(result[key], value)
            else:
                result[key] = value
        return result

    def get_genomic_targeting_config(self) -> dict[str, Any]:
        """
        Get genomic targeting configuration section.

        Returns:
            Genomic targeting configuration dictionary
        """
        return self.get_nested("genomic_targeting", default={})

    def is_genomic_targeting_enabled(self) -> bool:
        """
        Check if genomic targeting flags are enabled.

        Returns:
            True if genomic targeting is enabled
        """
        return self.get_nested("genomic_targeting", "enabled", default=False)

    def get_genomic_targeting_file_path(self) -> str:
        """
        Get genomic targeting file path.

        Returns:
            File path string for genomic targeting flags
        """
        return self.get_nested(
            "genomic_targeting",
            "file_path",
            default="data/manual/genomic_targeting_flags.xlsx",
        )

    def get_genomic_targeting_gene_column(self) -> str:
        """
        Get column name for gene symbols in targeting file.

        Returns:
            Column name string for gene symbols
        """
        return self.get_nested(
            "genomic_targeting", "gene_column", default="gene_symbol"
        )

    def get_genomic_targeting_column(self) -> str:
        """
        Get column name for targeting flags in targeting file.

        Returns:
            Column name string for targeting flags
        """
        return self.get_nested(
            "genomic_targeting", "targeting_column", default="targeting"
        )

    def get_genomic_targeting_default_value(self) -> bool:
        """
        Get default value for genes not found in targeting file.

        Returns:
            Default targeting flag value (boolean)
        """
        return self.get_nested("genomic_targeting", "default_value", default=False)

    def is_genomic_targeting_missing_file_allowed(self) -> bool:
        """
        Check if missing targeting file is allowed.

        Returns:
            True if missing targeting file should not cause failure
        """
        return self.get_nested("genomic_targeting", "allow_missing_file", default=True)

    def is_genomic_targeting_validation_enabled(self) -> bool:
        """
        Check if gene symbol validation is enabled for targeting file.

        Returns:
            True if gene symbols should be validated against HGNC
        """
        return self.get_nested(
            "genomic_targeting", "validate_gene_symbols", default=False
        )

    def to_dict(self) -> dict[str, Any]:
        """
        Get the full configuration as a dictionary.

        Returns:
            Complete configuration dictionary
        """
        return self.config.copy()
