"""
Centogene parser.

This module extracts gene information from Centogene panel pages.
Specifically targets genes in modal tables.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class CentogeneParser(BaseParser):
    """Parser for Centogene panels."""

    def parse(self) -> list[str]:
        """
        Parse Centogene panel page and extract gene symbols.

        Returns:
            List of gene symbols found on the page

        Raises:
            Exception: If parsing fails
        """
        try:
            # Make request with specific network error handling
            try:
                response = requests.get(self.url, timeout=30)
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(
                    f"Network request to Centogene URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Find modal with ID containing "74420" and look for table
            modal = soup.find("div", id=lambda x: x and "74420" in x)
            if modal:
                table = modal.find("table")
                if table:
                    # Extract genes from table cells
                    cells = table.find_all(["td", "th"])
                    for cell in cells:
                        text = cell.get_text(strip=True)
                        if text and len(text) <= 20:
                            cleaned_gene = self.clean_gene_symbol(text)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # If modal approach didn't work, look for tables generally
            if not genes:
                tables = soup.find_all("table")
                for table in tables:
                    cells = table.find_all(["td", "th"])
                    for cell in cells:
                        text = cell.get_text(strip=True)
                        if text and len(text) <= 20:
                            cleaned_gene = self.clean_gene_symbol(text)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # Look for gene information in other elements
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower() for keyword in ["gene", "panel", "test"]
                    ):
                        import re

                        potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,7}\b", text)
                        for gene in potential_genes:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # Remove duplicates while preserving order
            unique_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)

            logger.info(f"Extracted {len(unique_genes)} genes from Centogene panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Centogene panel: {e}")
            raise
