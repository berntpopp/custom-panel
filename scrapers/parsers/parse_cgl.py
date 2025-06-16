"""
CGL (Cincinnati Children's) parser.

This module extracts gene information from CGL panel pages.
Looks for genes in paragraphs following "GENES TARGETED" headers.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class CglParser(BaseParser):
    """Parser for CGL (Cincinnati Children's) panels."""

    def parse(self) -> list[str]:
        """
        Parse CGL panel page and extract gene symbols.

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

            # Find h2 with "GENES TARGETED" text and look for following p tags
            genes_targeted_headers = soup.find_all(
                "h2", string=lambda text: text and "genes targeted" in text.lower()
            )

            for header in genes_targeted_headers:
                # Find the next paragraph
                current = header.next_sibling
                while current:
                    if hasattr(current, "name") and current.name == "p":
                        text = current.get_text(strip=True)
                        if text:
                            # Parse comma-separated gene list with significant string manipulation
                            # Split by various delimiters
                            for delimiter in [",", ";", "\n", "\t"]:
                                parts = text.split(delimiter)
                                for part in parts:
                                    cleaned_gene = self.clean_gene_symbol(part)
                                    if cleaned_gene and self.validate_gene_symbol(
                                        cleaned_gene
                                    ):
                                        genes.append(cleaned_gene)
                        break
                    current = current.next_sibling

            # If specific approach didn't work, look for gene patterns generally
            if not genes:
                for element in soup.find_all(["p", "div", "span", "li"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower()
                        for keyword in ["gene", "panel", "target"]
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

            logger.info(f"Extracted {len(unique_genes)} genes from CGL panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing CGL panel: {e}")
            raise
