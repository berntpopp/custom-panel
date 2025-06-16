"""
MGZ Munich parser.

This module extracts gene information from MGZ Munich panel pages.
Looks for genes in divs with class="panel_gene".
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class MgzParser(BaseParser):
    """Parser for MGZ Munich panels."""

    def parse(self) -> list[str]:
        """
        Parse MGZ panel page and extract gene symbols.

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

            genes = []

            # Find divs with class="panel_gene"
            panel_gene_divs = soup.find_all("div", class_="panel_gene")
            for div in panel_gene_divs:
                text = div.get_text(strip=True)
                if text:
                    cleaned_gene = self.clean_gene_symbol(text)
                    if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                        genes.append(cleaned_gene)

            # If no specific divs found, look for other gene indicators
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li"]):
                    class_attr = element.get("class", [])
                    if any("gene" in str(cls).lower() for cls in class_attr):
                        text = element.get_text(strip=True)
                        if text and len(text) <= 20:
                            cleaned_gene = self.clean_gene_symbol(text)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # General search for gene patterns
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower()
                        for keyword in ["gene", "panel", "diagnostic"]
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

            logger.info(f"Extracted {len(unique_genes)} genes from MGZ panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing MGZ panel: {e}")
            raise
