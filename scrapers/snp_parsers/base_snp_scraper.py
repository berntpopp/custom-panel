"""
Abstract base class for all SNP panel scrapers.

This module defines the BaseSNPScraper class that all individual SNP scrapers must inherit from.
"""

import logging
import re
from abc import ABC, abstractmethod
from typing import Any

import requests
from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)


class BaseSNPScraper(ABC):
    """
    Abstract base class for SNP panel scrapers.

    All SNP scraper implementations must inherit from this class and implement the parse method.
    """

    def __init__(self, url: str, config: dict[str, Any] | None = None) -> None:
        """
        Initialize the scraper with URL and optional configuration.

        Args:
            url: The URL to scrape data from
            config: Optional configuration dictionary
        """
        self.url = url
        self.config = config or {}
        self.session = requests.Session()
        self.session.headers.update(
            {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
            },
        )

    @abstractmethod
    def parse(self) -> dict[str, Any]:
        """
        Parse the target URL and extract SNP data.

        Returns:
            Dictionary containing:
                - panel_name: Name of the panel
                - source_url: Source URL
                - snps: List of dictionaries with SNP data
                - metadata: Optional metadata about the panel

        Raises:
            Exception: If parsing fails
        """

    def fetch_content(self, url: str = None, use_selenium: bool = False) -> str:
        """
        Fetch content from URL using requests or selenium.

        Args:
            url: URL to fetch (defaults to self.url)
            use_selenium: Whether to use selenium for JavaScript-rendered content

        Returns:
            Page content as string
        """
        url = url or self.url

        if use_selenium:
            # Import here to avoid dependency if not needed
            from custom_panel.sources_snp.downloader import PanelDownloader

            with PanelDownloader() as downloader:
                file_path = downloader.download(url, "html", "selenium")
                with open(file_path, encoding="utf-8") as f:
                    return f.read()
        else:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            return response.text

    def clean_rsid(self, rsid: str) -> str:
        """
        Clean and standardize an rsID.

        Args:
            rsid: Raw rsID string

        Returns:
            Cleaned rsID
        """
        if not rsid:
            return ""

        # Remove whitespace
        rsid = rsid.strip()

        # Ensure rs prefix is lowercase
        if rsid.upper().startswith("RS"):
            rsid = "rs" + rsid[2:]

        # Remove any non-alphanumeric characters after the rsID
        match = re.match(r"(rs\d+)", rsid)
        if match:
            return match.group(1)

        return rsid

    def validate_rsid(self, rsid: str) -> bool:
        """
        Validate that a string is a valid rsID.

        Args:
            rsid: rsID to validate

        Returns:
            True if the rsID appears valid
        """
        if not rsid:
            return False

        # Must match pattern rs followed by digits
        return bool(re.match(r"^rs\d+$", rsid))

    def extract_rsids_from_text(self, text: str) -> list[str]:
        """
        Extract all valid rsIDs from a text string.

        Args:
            text: Text to search for rsIDs

        Returns:
            List of unique, validated rsIDs
        """
        # Find all potential rsIDs
        pattern = r"\b(rs\d+)\b"
        matches = re.findall(pattern, text, re.IGNORECASE)

        # Clean and validate
        rsids = []
        seen = set()
        for match in matches:
            rsid = self.clean_rsid(match)
            if self.validate_rsid(rsid) and rsid not in seen:
                seen.add(rsid)
                rsids.append(rsid)

        return rsids

    def parse_table_for_rsids(
        self,
        soup: BeautifulSoup,
        table_selector: str = None,
        rsid_column: str = None,
    ) -> list[str]:
        """
        Extract rsIDs from HTML tables.

        Args:
            soup: BeautifulSoup object
            table_selector: CSS selector for tables
            rsid_column: Column name containing rsIDs

        Returns:
            List of rsIDs found in tables
        """
        rsids = []

        # Find tables
        if table_selector:
            tables = soup.select(table_selector)
        else:
            tables = soup.find_all("table")

        for table in tables:
            try:
                # Convert to pandas DataFrame
                from io import StringIO

                import pandas as pd

                df_list = pd.read_html(StringIO(str(table)))
                if not df_list:
                    continue

                df = df_list[0]

                # Look for rsID column
                if rsid_column and rsid_column in df.columns:
                    col_data = df[rsid_column].dropna().astype(str)
                    for value in col_data:
                        rsids.extend(self.extract_rsids_from_text(value))
                else:
                    # Search all columns
                    for col in df.columns:
                        col_data = df[col].dropna().astype(str)
                        for value in col_data:
                            rsids.extend(self.extract_rsids_from_text(value))

            except Exception as e:
                logger.warning(f"Failed to parse table: {e}")
                continue

        return list(set(rsids))  # Return unique rsIDs

    def create_snp_record(
        self,
        rsid: str,
        chromosome: str = None,
        position: int = None,
        start: int = None,
        end: int = None,
        strand: str = None,
        assembly: str = None,
        **kwargs,
    ) -> dict[str, Any]:
        """
        Create a standardized SNP record with optional genomic coordinates.

        Args:
            rsid: The rsID
            chromosome: Chromosome (e.g., '1', 'chr1', 'X')
            position: Genomic position (for point mutations)
            start: Start position (for ranges)
            end: End position (for ranges)
            strand: Strand orientation ('+', '-')
            assembly: Reference assembly ('GRCh37', 'GRCh38', 'hg19', 'hg38')
            **kwargs: Additional fields for the record

        Returns:
            Standardized SNP record dictionary
        """
        record = {"rsid": rsid, "source": self.config.get("name", "unknown")}

        # Add genomic coordinates if provided
        if chromosome is not None:
            record["chromosome"] = self.clean_chromosome(chromosome)

        if position is not None:
            record["position"] = int(position)
        elif start is not None:
            record["start"] = int(start)
            if end is not None:
                record["end"] = int(end)

        if strand is not None:
            record["strand"] = strand

        if assembly is not None:
            record["assembly"] = self.clean_assembly(assembly)

        record.update(kwargs)
        return record

    def clean_chromosome(self, chromosome: str) -> str:
        """
        Clean and standardize chromosome notation.

        Args:
            chromosome: Raw chromosome string

        Returns:
            Cleaned chromosome (e.g., '1', '2', 'X', 'Y', 'MT')
        """
        if not chromosome:
            return ""

        chrom = str(chromosome).strip().upper()

        # Remove 'chr' prefix if present
        if chrom.startswith("CHR"):
            chrom = chrom[3:]

        # Handle mitochondrial
        if chrom in ["MT", "M", "MITO", "MITOCHONDRIAL"]:
            return "MT"

        # Handle sex chromosomes
        if chrom in ["23", "XX"]:
            return "X"
        if chrom in ["24", "XY"]:
            return "Y"

        # Numeric chromosomes
        if chrom.isdigit():
            return chrom

        # Already clean format
        if chrom in ["X", "Y"]:
            return chrom

        return chrom

    def clean_assembly(self, assembly: str) -> str:
        """
        Clean and standardize genome assembly notation.

        Args:
            assembly: Raw assembly string

        Returns:
            Standardized assembly ('GRCh37', 'GRCh38')
        """
        if not assembly:
            return ""

        assembly = str(assembly).strip().lower()

        # Map common assembly names
        assembly_map = {
            "grch37": "GRCh37",
            "grch38": "GRCh38",
            "hg19": "GRCh37",
            "hg38": "GRCh38",
            "b37": "GRCh37",
            "b38": "GRCh38",
            "37": "GRCh37",
            "38": "GRCh38",
        }

        return assembly_map.get(assembly, assembly.upper())

    def parse_genomic_position(self, position_str: str) -> dict[str, Any]:
        """
        Parse genomic position string in various formats.

        Args:
            position_str: Position string (e.g., 'chr1:123456', '1:123456-123457')

        Returns:
            Dictionary with parsed coordinates
        """
        result = {}

        if not position_str:
            return result

        position_str = str(position_str).strip()

        # Pattern: chr1:123456 or chr1:123456-123457
        match = re.match(r"(chr)?(\w+):(\d+)(?:-(\d+))?", position_str, re.IGNORECASE)
        if match:
            chromosome = match.group(2)
            start = int(match.group(3))
            end = int(match.group(4)) if match.group(4) else start

            result["chromosome"] = self.clean_chromosome(chromosome)
            result["start"] = start
            result["end"] = end

            return result

        # Pattern: just a number (position only)
        if position_str.isdigit():
            result["position"] = int(position_str)

        return result
