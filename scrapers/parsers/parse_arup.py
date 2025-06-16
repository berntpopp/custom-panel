"""
ARUP Laboratories parser.

This module extracts gene information from ARUP panel pages.
Looks for genes in tables following elements with id="genes-tested1".
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class ArupParser(BaseParser):
    """Parser for ARUP Laboratories panels."""

    def parse(self) -> list[str]:
        """
        Parse ARUP panel page and extract gene symbols.

        Returns:
            List of gene symbols found on the page

        Raises:
            Exception: If parsing fails
        """
        try:
            # Make request
            try:
                response = requests.get(self.url, timeout=30)
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(f"Network request to ARUP URL failed: {self.url} - {e}")
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Find element with id="genes-tested1" and look for following table
            genes_tested_element = soup.find(id="genes-tested1")
            if genes_tested_element:
                # Find the next table
                current = genes_tested_element.next_sibling
                while current:
                    if hasattr(current, "name") and current.name == "table":
                        # Extract genes from this table
                        cells = current.find_all(["td", "th"])
                        for cell in cells:
                            text = cell.get_text(strip=True)
                            if text and len(text) <= 20:
                                cleaned_gene = self.clean_gene_symbol(text)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene
                                ):
                                    genes.append(cleaned_gene)
                        break
                    current = current.next_sibling

            # If specific approach didn't work, look for genes in all tables
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

            # General search for gene patterns
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

            logger.info(f"Extracted {len(unique_genes)} genes from ARUP panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing ARUP panel: {e}")
            raise
