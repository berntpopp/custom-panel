"""
University of Chicago parser.

This module extracts gene information from University of Chicago panel pages.
Looks for genes in links with href containing "/gene/".
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class UChicagoParser(BaseParser):
    """Parser for University of Chicago panels."""

    def parse(self) -> list[str]:
        """
        Parse University of Chicago panel page and extract gene symbols.

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
                logger.error(
                    f"Network request to University of Chicago URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Find links with href containing "/gene/"
            gene_links = soup.find_all("a", href=lambda x: x and "/gene/" in x)
            for link in gene_links:
                text = link.get_text(strip=True)
                if text:
                    cleaned_gene = self.clean_gene_symbol(text)
                    if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                        genes.append(cleaned_gene)

            # If no gene links found, look for other patterns
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li", "td"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower()
                        for keyword in ["gene", "panel", "cancer"]
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

            logger.info(
                f"Extracted {len(unique_genes)} genes from University of Chicago panel"
            )
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing University of Chicago panel: {e}")
            raise
