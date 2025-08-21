"""
Output manager for intermediate files and structured logging.

This module provides functionality to save intermediate files at each pipeline step
and manage structured logging for debugging purposes.
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from .io import save_panel_data

logger = logging.getLogger(__name__)


class OutputManager:
    """Manages intermediate file output and structured logging."""

    def __init__(self, config: dict[str, Any], output_dir: str | Path):
        """
        Initialize the output manager.

        Args:
            config: Configuration dictionary
            output_dir: Base output directory
        """
        self.config = config
        self.base_output_dir = Path(output_dir)

        # Get configuration settings
        self.output_config = config.get("output", {})
        self.intermediate_config = self.output_config.get("intermediate_files", {})
        self.logging_config = self.output_config.get("file_logging", {})
        self.dir_config = config.get("directory_structure", {})

        # Check if intermediate files are enabled
        self.intermediate_enabled = self.intermediate_config.get("enabled", False)
        self.file_logging_enabled = self.logging_config.get("enabled", False)

        # Debug: Log the configuration
        logger.info(
            f"OutputManager initialized with intermediate_enabled={self.intermediate_enabled}, format={self.intermediate_config.get('format', 'NOT_SET')}",
        )
        logger.debug(f"Full intermediate config: {self.intermediate_config}")

        # Set up directory structure
        self.use_structured = self.dir_config.get("use_structured_output", True)
        self.subdirs = self.dir_config.get("subdirs", {})
        self.timestamp_format = self.dir_config.get("timestamp_format", "%Y%m%d_%H%M%S")

        # Create run-specific directory if using structured output
        if self.use_structured:
            timestamp = datetime.now().strftime(self.timestamp_format)
            self.run_dir = self.base_output_dir / f"run_{timestamp}"
        else:
            self.run_dir = self.base_output_dir

        # Initialize logging if enabled
        if self.file_logging_enabled:
            self._setup_file_logging()

        logger.info(f"Output manager initialized - Run directory: {self.run_dir}")

    def _setup_file_logging(self) -> None:
        """Set up file-based logging."""
        log_dir = self.run_dir / self.subdirs.get("logs", "logs")
        log_dir.mkdir(parents=True, exist_ok=True)

        log_format = self.logging_config.get("format", "text")
        include_debug = self.logging_config.get("include_debug", True)
        separate_components = self.logging_config.get("separate_by_component", True)

        # Set log level
        log_level = logging.DEBUG if include_debug else logging.INFO

        if separate_components:
            # Create separate loggers for each component
            components = self.logging_config.get("components", [])
            for component in components:
                self._setup_component_logger(component, log_dir, log_format, log_level)

        # Always create a main pipeline log
        self._setup_component_logger("pipeline", log_dir, log_format, log_level)

    def _setup_component_logger(
        self, component: str, log_dir: Path, log_format: str, log_level: int,
    ) -> None:
        """Set up logging for a specific component."""
        log_file = log_dir / f"{component}.log"

        # Create file handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)

        # Set format
        if log_format == "json":
            formatter: logging.Formatter = JsonFormatter()
        else:
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            )
        file_handler.setFormatter(formatter)

        # Add handler to the specific logger
        component_logger = logging.getLogger(f"custom_panel.{component}")
        component_logger.addHandler(file_handler)

        logger.debug(f"Set up file logging for {component}: {log_file}")

    def save_intermediate_data(
        self,
        data: pd.DataFrame,
        step: str,
        description: str = "",
        metadata: dict[str, Any] | None = None,
    ) -> Path | None:
        """
        Save intermediate data to file.

        Args:
            data: DataFrame to save
            step: Pipeline step name
            description: Human-readable description
            metadata: Additional metadata to save

        Returns:
            Path to saved file, or None if intermediate files are disabled
        """
        if not self.intermediate_enabled or data.empty:
            return None

        # Determine subdirectory based on step
        subdir_mapping = {
            "raw_data": self.subdirs.get("raw_data", "01_raw_data"),
            "standardized_data": self.subdirs.get(
                "standardized_data", "02_standardized_data",
            ),
            "merged_data": self.subdirs.get("merged_data", "03_merged_data"),
            "scored_data": self.subdirs.get("scored_data", "04_scored_data"),
            "annotated_data": self.subdirs.get("annotated_data", "05_annotated_data"),
        }

        step_dir = self.run_dir / subdir_mapping.get(step, step)
        step_dir.mkdir(parents=True, exist_ok=True)

        # Generate filename
        timestamp = datetime.now().strftime("%H%M%S")
        safe_description = (
            description.replace(" ", "_").replace("/", "_") if description else ""
        )
        if safe_description:
            filename = f"{timestamp}_{step}_{safe_description}"
        else:
            filename = f"{timestamp}_{step}"

        # Get file format
        file_format = self.intermediate_config.get("format", "csv")
        compress = self.intermediate_config.get("compress", False)

        # Debug: Log the format being used
        logger.debug(f"Using file_format='{file_format}' for step='{step}'")

        # Add compression extension if needed
        if compress and file_format in ["csv", "excel"]:
            if file_format == "csv":
                file_extension = "csv.gz"
            else:
                file_extension = "xlsx"  # Excel files are already compressed
        else:
            file_extension = file_format.replace("excel", "xlsx")

        file_path = step_dir / f"{filename}.{file_extension}"

        try:
            # Save the data - skip validation for data with different schema
            if step in ["scored_data", "annotated_data"] or (
                metadata and metadata.get("data_type") == "snp"
            ):
                # Scored, annotated, and SNP data have different schema, save directly
                self._save_data_direct(data, file_path, file_format)
            else:
                # Use standard validation for other data types
                save_panel_data(data, file_path, file_format)

            # Save metadata if provided
            if metadata:
                metadata_path = step_dir / f"{filename}_metadata.json"
                with open(metadata_path, "w") as f:
                    json.dump(metadata, f, indent=2, default=str)

            logger.info(f"Saved intermediate data: {file_path} ({len(data)} records)")
            return file_path

        except Exception as e:
            logger.error(f"Failed to save intermediate data to {file_path}: {e}")
            return None

    def save_source_data(self, data: pd.DataFrame, source_name: str) -> Path | None:
        """Save raw source data."""
        if not self.intermediate_config.get("include_raw_data", True):
            return None

        return self.save_intermediate_data(
            data, "raw_data", f"source_{source_name}", {"source": source_name},
        )

    def save_snp_data(self, data: pd.DataFrame, category: str) -> Path | None:
        """Save SNP data with appropriate handling for different schema."""
        if not self.intermediate_config.get("include_raw_data", True):
            return None

        # SNP data has a different schema than gene data, so we bypass validation
        return self.save_intermediate_data(
            data,
            "raw_data",
            f"SNP_{category}",
            {"category": category, "data_type": "snp"},
        )

    def save_regions_data(self, data: pd.DataFrame, region_type: str) -> Path | None:
        """Save regions data with appropriate handling for different schema."""
        if not self.intermediate_config.get("include_raw_data", True):
            return None

        # Regions data has a different schema than gene data, so we bypass validation
        return self.save_intermediate_data(
            data,
            "raw_data",
            f"REGIONS_{region_type}",
            {"region_type": region_type, "data_type": "regions"},
        )

    def save_standardized_data(
        self,
        data: pd.DataFrame,
        source_name: str,
        symbol_changes: dict[str, dict[str, str | None]],
    ) -> Path | None:
        """Save standardized source data."""
        if not self.intermediate_config.get("include_standardized_data", True):
            return None

        # Count actual changes
        changes_count = sum(
            1
            for k, v in symbol_changes.items()
            if v["approved_symbol"] is not None and k != v["approved_symbol"]
        )

        metadata = {
            "source": source_name,
            "symbol_changes": symbol_changes,
            "changes_count": changes_count,
        }

        return self.save_intermediate_data(
            data, "standardized_data", f"standardized_{source_name}", metadata,
        )

    def save_merged_data(
        self, data: pd.DataFrame, source_stats: dict[str, int],
    ) -> Path | None:
        """Save merged data before scoring."""
        if not self.intermediate_config.get("include_merged_data", True):
            return None

        metadata = {
            "source_statistics": source_stats,
            "total_records": len(data),
            "unique_genes": (
                data["approved_symbol"].nunique()
                if "approved_symbol" in data.columns
                else 0
            ),
        }

        return self.save_intermediate_data(
            data, "merged_data", "merged_all_sources", metadata,
        )

    def save_scored_data(
        self, data: pd.DataFrame, scoring_summary: dict[str, Any],
    ) -> Path | None:
        """Save scored data before decision logic."""
        if not self.intermediate_config.get("include_scored_data", True):
            return None

        return self.save_intermediate_data(
            data, "scored_data", "scored_genes", scoring_summary,
        )

    def save_annotated_data(
        self, data: pd.DataFrame, annotation_summary: dict[str, Any],
    ) -> Path | None:
        """Save final annotated data."""
        if not self.intermediate_config.get("include_annotated_data", True):
            return None

        return self.save_intermediate_data(
            data, "annotated_data", "final_annotated", annotation_summary,
        )

    def get_final_output_dir(self) -> Path:
        """Get the directory for final output files."""
        final_dir = self.run_dir / self.subdirs.get("final_output", "06_final_output")
        final_dir.mkdir(parents=True, exist_ok=True)
        return final_dir

    def get_run_summary(self) -> dict[str, Any]:
        """Get summary of the current run."""
        return {
            "run_directory": str(self.run_dir),
            "timestamp": datetime.now().isoformat(),
            "intermediate_files_enabled": self.intermediate_enabled,
            "file_logging_enabled": self.file_logging_enabled,
            "output_format": self.intermediate_config.get("format", "csv"),
            "structured_output": self.use_structured,
        }

    def cleanup_old_runs(self, keep_runs: int = 10) -> None:
        """Clean up old run directories to save space."""
        if not self.use_structured:
            return

        try:
            # Find all run directories
            run_dirs = [
                d
                for d in self.base_output_dir.iterdir()
                if d.is_dir() and d.name.startswith("run_")
            ]

            # Sort by creation time (newest first)
            run_dirs.sort(key=lambda x: x.stat().st_ctime, reverse=True)

            # Remove old runs beyond the keep limit
            for old_run in run_dirs[keep_runs:]:
                try:
                    import shutil

                    shutil.rmtree(old_run)
                    logger.info(f"Cleaned up old run directory: {old_run}")
                except Exception as e:
                    logger.warning(f"Failed to clean up {old_run}: {e}")

        except Exception as e:
            logger.warning(f"Failed to cleanup old runs: {e}")

    def _save_data_direct(self, df: pd.DataFrame, path: Path, format: str) -> None:
        """
        Save data directly without validation (for scored data with different schema).

        Args:
            df: DataFrame to save
            path: Output file path
            format: Output format
        """
        path.parent.mkdir(parents=True, exist_ok=True)

        if format.lower() == "parquet":
            df.to_parquet(path, index=False, engine="pyarrow")
        elif format.lower() == "csv":
            df.to_csv(path, index=False)
        elif format.lower() == "excel":
            df.to_excel(path, index=False, engine="openpyxl")
        else:
            raise ValueError(f"Unsupported format: {format}")

        logger.info(f"Saved {len(df)} records to {path} (direct save)")


class JsonFormatter(logging.Formatter):
    """JSON formatter for structured logging."""

    def format(self, record: logging.LogRecord) -> str:
        log_obj = {
            "timestamp": datetime.fromtimestamp(record.created).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }

        if record.exc_info:
            log_obj["exception"] = self.formatException(record.exc_info)

        return json.dumps(log_obj)
