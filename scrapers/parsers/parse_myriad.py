"""
Myriad Genetics myRisk parser.

This module extracts gene information from Myriad Genetics myRisk panel pages.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class MyriadParser(BaseParser):
    """Parser for Myriad Genetics myRisk panels."""

    def parse(self) -> list[str]:
        """
        Parse Myriad myRisk page and extract gene symbols.

        Returns:
            List of gene symbols found on the page

        Raises:
            Exception: If parsing fails
        """
        try:
            # Make request
            response = requests.get(self.url, timeout=30)
            response.raise_for_status()

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            # Find gene cells using XPath equivalent in BeautifulSoup
            # Target table cells with class="gene"
            gene_cells = soup.find_all("td", class_=lambda x: x and "gene" in x.lower())

            genes = []
            for cell in gene_cells:
                gene_text = cell.get_text(strip=True)
                if gene_text:
                    # Clean gene symbol (remove parenthetical content)
                    cleaned_gene = self.clean_gene_symbol(gene_text)
                    if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                        genes.append(cleaned_gene)

            # Also try to find genes in other common patterns
            if not genes:
                # Look for gene lists in paragraphs or divs
                for element in soup.find_all(["p", "div", "span"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower() for keyword in ["gene", "panel", "test"]
                    ):
                        # Extract potential gene symbols (uppercase 2-8 character words)
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

            logger.info(f"Extracted {len(unique_genes)} genes from Myriad myRisk panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Myriad myRisk panel: {e}")
            raise
