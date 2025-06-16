"""
Fulgent Genetics parser.

This module extracts gene information from Fulgent Genetics panel pages.
Specifically looks for gene lists in meta description tags.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class FulgentParser(BaseParser):
    """Parser for Fulgent Genetics panels."""

    def parse(self) -> list[str]:
        """
        Parse Fulgent Genetics page and extract gene symbols.

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
                    f"Network request to Fulgent Genetics URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # R script uses: //meta[contains(@content,"AIP")]
            # Find meta tags that contain "AIP" in their content
            import re

            meta_tags = soup.find_all("meta")
            for meta in meta_tags:
                content = meta.get("content", "")
                if "AIP" in content:
                    # Parse genes from content using R script approach

                    # Clean and split the content
                    content = (
                        content.replace("Full", "")
                        .replace("Comprehensive", "")
                        .replace("Cancer", "")
                        .replace("Panel", "")
                        .replace("COMPREHENSIVE", "")
                    )

                    # Split by common separators and extract gene symbols
                    potential_genes = re.split(r"[,;\s]+", content)
                    for gene in potential_genes:
                        gene = gene.strip()
                        if gene:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                    break

            # If no specific pattern found, extract all potential gene symbols
            if not genes:
                for meta in meta_tags:
                    content = meta.get("content", "")
                    if content:
                        potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,7}\b", content)
                        for gene in potential_genes:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # If meta description didn't work, look in page content
            if not genes:
                # Look for gene information in the page body
                for element in soup.find_all(["div", "p", "span", "td", "li"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower()
                        for keyword in ["gene", "panel", "target"]
                    ):
                        # Extract potential gene symbols
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
                f"Extracted {len(unique_genes)} genes from Fulgent Genetics panel"
            )
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Fulgent Genetics panel: {e}")
            raise
