"""
Mayo Clinic parser.

This module extracts gene information from Mayo Clinic panel pages.
Looks for genes in italic tags.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class MayoParser(BaseParser):
    """Parser for Mayo Clinic panels."""

    def parse(self) -> list[str]:
        """
        Parse Mayo Clinic panel page and extract gene symbols.

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
                    f"Network request to Mayo Clinic URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Find italic tags containing gene lists
            italic_tags = soup.find_all("i")
            for italic in italic_tags:
                text = italic.get_text(strip=True)
                if text and (
                    "," in text or len(text.split()) > 1
                ):  # Likely a gene list
                    # Split by comma and extract genes
                    for gene_part in text.split(","):
                        cleaned_gene = self.clean_gene_symbol(gene_part)
                        if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                            genes.append(cleaned_gene)

            # Also check em tags
            if not genes:
                em_tags = soup.find_all("em")
                for em in em_tags:
                    text = em.get_text(strip=True)
                    if text and ("," in text or len(text.split()) > 1):
                        for gene_part in text.split(","):
                            cleaned_gene = self.clean_gene_symbol(gene_part)
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

            logger.info(f"Extracted {len(unique_genes)} genes from Mayo Clinic panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Mayo Clinic panel: {e}")
            raise
