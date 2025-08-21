#!/usr/bin/env python3
"""
Master runner script for commercial panel scrapers.

This script coordinates the execution of individual parser classes to fetch
gene data from commercial diagnostic panel websites and save the results
in a standardized JSON format.
"""

import argparse
import importlib
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

logger = logging.getLogger(__name__)


def load_config(config_path: str | None = None) -> dict[str, Any]:
    """
    Load configuration from YAML file.

    Args:
        config_path: Path to configuration file. If None, uses default.

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config file is invalid
    """
    if config_path is None:
        # Use default config
        default_config = (
            Path(__file__).parent.parent
            / "custom_panel"
            / "config"
            / "default_config.yml"
        )
        config_path = str(default_config)

    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_file) as f:
        config = yaml.safe_load(f)

    return config


def get_parser_class(parser_module: str, parser_class: str) -> type:
    """
    Dynamically import and return a parser class.

    Args:
        parser_module: Module name (e.g., 'parse_myriad')
        parser_class: Class name (e.g., 'MyriadParser')

    Returns:
        Parser class

    Raises:
        ImportError: If module or class cannot be imported
    """
    try:
        module = importlib.import_module(f"scrapers.parsers.{parser_module}")
        parser_cls = getattr(module, parser_class)

        # Add parent directory to path to import custom_panel modules
        sys.path.append(str(Path(__file__).parent.parent))
        from scrapers.parsers.base_parser import BaseParser

        if not issubclass(parser_cls, BaseParser):
            raise TypeError(f"{parser_class} must inherit from BaseParser")

        return parser_cls

    except (ImportError, AttributeError) as e:
        raise ImportError(
            f"Failed to import {parser_class} from {parser_module}: {e}",
        ) from e


def create_output_json(
    panel_name: str,
    source_url: str,
    genes: list[str],
    output_path: str | Path,
) -> None:
    """
    Create standardized JSON output file.

    Args:
        panel_name: Name of the panel
        source_url: URL the data was retrieved from
        genes: List of gene symbols
        output_path: Path where JSON file should be saved
    """
    output_data = {
        "panel_name": panel_name,
        "source_url": source_url,
        "retrieval_date": datetime.now().strftime("%Y-%m-%d"),
        "genes": sorted(set(genes)),  # Remove duplicates and sort
    }

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        json.dump(output_data, f, indent=2)

    logger.info(f"Saved {len(output_data['genes'])} genes to {output_file}")


def run_scraper(scraper_name: str, scraper_config: dict[str, Any]) -> None:
    """
    Run a single scraper and save its output.

    Args:
        scraper_name: Name of the scraper
        scraper_config: Configuration for this scraper

    Raises:
        Exception: If scraping fails
    """
    logger.info(f"Running scraper: {scraper_name}")

    try:
        # Get configuration
        url = scraper_config["url"]
        parser_module = scraper_config["parser_module"]
        parser_class = scraper_config["parser_class"]
        output_path = scraper_config["output_path"]

        # Handle multiple URLs for some scrapers (like Blueprint)
        urls = [url] if isinstance(url, str) else url

        # Handle subpanel URLs for Blueprint Genetics
        if "subpanel_urls" in scraper_config:
            urls.extend(scraper_config["subpanel_urls"])

        all_genes = []

        for single_url in urls:
            # Get parser class and instantiate
            parser_cls = get_parser_class(parser_module, parser_class)
            parser = parser_cls(single_url, scraper_config)

            # Parse genes
            genes = parser.parse()
            all_genes.extend(genes)

            logger.info(f"Extracted {len(genes)} genes from {single_url}")

        # Create output JSON
        create_output_json(
            panel_name=scraper_name,
            source_url=urls[0] if len(urls) == 1 else str(urls),
            genes=all_genes,
            output_path=output_path,
        )

        logger.info(f"Successfully completed scraper: {scraper_name}")

    except Exception as e:
        logger.error(f"Failed to run scraper {scraper_name}: {e}")
        raise


def main() -> None:
    """Main entry point for the scraper runner."""
    parser = argparse.ArgumentParser(
        description="Run commercial panel scrapers to extract gene data",
    )
    parser.add_argument(
        "--config",
        type=str,
        help="Path to configuration file (default: uses default config)",
    )
    parser.add_argument(
        "--names",
        type=str,
        nargs="+",
        help="Specific scraper names to run (default: run all enabled scrapers)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        help="Override output directory for all scrapers",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without actually running scrapers",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        # Load configuration
        config = load_config(args.config)
        scrapers_config = config.get("scrapers", {})

        if not scrapers_config:
            logger.error("No scrapers configuration found in config file")
            sys.exit(1)

        # Determine which scrapers to run
        if args.names:
            scrapers_to_run = {
                name: scrapers_config[name]
                for name in args.names
                if name in scrapers_config
            }
            missing = set(args.names) - set(scrapers_config.keys())
            if missing:
                logger.warning(f"Unknown scrapers: {', '.join(missing)}")
        else:
            scrapers_to_run = {
                name: config
                for name, config in scrapers_config.items()
                if config.get("enabled", True)
            }

        if not scrapers_to_run:
            logger.error("No scrapers to run")
            sys.exit(1)

        logger.info(
            f"Will run {len(scrapers_to_run)} scrapers: {', '.join(scrapers_to_run.keys())}",
        )

        if args.dry_run:
            logger.info("Dry run mode - no actual scraping will be performed")
            for name, config in scrapers_to_run.items():
                output_path = args.output_dir or config["output_path"]
                logger.info(f"  {name}: {config['url']} -> {output_path}")
            return

        # Run scrapers
        success_count = 0
        for name, scraper_config in scrapers_to_run.items():
            try:
                # Override output directory if specified
                if args.output_dir:
                    output_file = Path(args.output_dir) / f"{name}.json"
                    scraper_config = scraper_config.copy()
                    scraper_config["output_path"] = str(output_file)

                run_scraper(name, scraper_config)
                success_count += 1

            except Exception as e:
                logger.error(f"Scraper {name} failed: {e}")
                continue

        logger.info(
            f"Completed {success_count}/{len(scrapers_to_run)} scrapers successfully",
        )

        if success_count == 0:
            sys.exit(1)

    except Exception as e:
        logger.error(f"Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
