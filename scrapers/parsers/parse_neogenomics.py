"""
NeoGenomics parser.

This module extracts gene information from NeoGenomics panel pages.
Looks for genes in paragraphs containing "AIP" and parses comma-separated lists.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class NeogenomicsParser(BaseParser):
    """Parser for NeoGenomics panels."""

    def parse(self) -> list[str]:
        """
        Parse NeoGenomics panel page and extract gene symbols.

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
                    f"Network request to NeoGenomics URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Find paragraphs containing "AIP"
            aip_paragraphs = soup.find_all(
                "p", string=lambda text: text and "AIP" in text
            )

            for paragraph in aip_paragraphs:
                text = paragraph.get_text(strip=True)
                # Parse comma-separated gene list
                for gene_part in text.split(","):
                    cleaned_gene = self.clean_gene_symbol(gene_part)
                    if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                        genes.append(cleaned_gene)

            # Also look for AIP in other elements
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li"]):
                    text = element.get_text()
                    if "AIP" in text:
                        for gene_part in text.split(","):
                            cleaned_gene = self.clean_gene_symbol(gene_part)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # General search for gene patterns if AIP approach didn't work
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li", "td"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower()
                        for keyword in ["gene", "panel", "tumor"]
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

            logger.info(f"Extracted {len(unique_genes)} genes from NeoGenomics panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing NeoGenomics panel: {e}")
            raise
