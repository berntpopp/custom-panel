"""
Abstract base class for all SNP panel scrapers.

This module defines the BaseSNPScraper class that all individual SNP scrapers must inherit from.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, List
import re
import logging
from pathlib import Path
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
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })

    @abstractmethod
    def parse(self) -> Dict[str, Any]:
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
        pass

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
                with open(file_path, 'r', encoding='utf-8') as f:
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
        match = re.match(r'(rs\d+)', rsid)
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
        return bool(re.match(r'^rs\d+$', rsid))

    def extract_rsids_from_text(self, text: str) -> List[str]:
        """
        Extract all valid rsIDs from a text string.

        Args:
            text: Text to search for rsIDs

        Returns:
            List of unique, validated rsIDs
        """
        # Find all potential rsIDs
        pattern = r'\b(rs\d+)\b'
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

    def parse_table_for_rsids(self, soup: BeautifulSoup, 
                             table_selector: str = None,
                             rsid_column: str = None) -> List[str]:
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
            tables = soup.find_all('table')
        
        for table in tables:
            try:
                # Convert to pandas DataFrame
                import pandas as pd
                from io import StringIO
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

    def create_snp_record(self, rsid: str, **kwargs) -> Dict[str, Any]:
        """
        Create a standardized SNP record.

        Args:
            rsid: The rsID
            **kwargs: Additional fields for the record

        Returns:
            Standardized SNP record dictionary
        """
        record = {
            "rsid": rsid,
            "source": self.config.get("name", "unknown")
        }
        record.update(kwargs)
        return record