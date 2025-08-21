"""
CEGAT parser.

This module extracts gene information from CEGAT panel pages.
Looks for genes in emphasis tags following Gene Directory headers.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class CegatParser(BaseParser):
    """Parser for CEGAT panels."""

    def parse(self) -> list[str]:
        """
        Parse CEGAT panel page and extract gene symbols.

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
                logger.error(f"Network request to CEGAT URL failed: {self.url} - {e}")
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Find h2 with "Gene Directory" text and look for following em tags
            gene_directory_headers = soup.find_all(
                "h2", string=lambda text: text and "gene directory" in text.lower(),
            )

            for header in gene_directory_headers:
                # Find all em tags that follow this header
                current = header.next_sibling
                while current:
                    if hasattr(current, "find_all"):
                        em_tags = current.find_all("em")
                        for em in em_tags:
                            text = em.get_text(strip=True)
                            if text:
                                # Check if this is a comma-separated list of genes
                                if "," in text:
                                    # Split by comma and process each gene
                                    potential_genes = [
                                        g.strip() for g in text.split(",")
                                    ]
                                    for gene in potential_genes:
                                        cleaned_gene = self.clean_gene_symbol(gene)
                                        if cleaned_gene and self.validate_gene_symbol(
                                            cleaned_gene,
                                        ):
                                            genes.append(cleaned_gene)
                                else:
                                    # Single gene
                                    cleaned_gene = self.clean_gene_symbol(text)
                                    if cleaned_gene and self.validate_gene_symbol(
                                        cleaned_gene,
                                    ):
                                        genes.append(cleaned_gene)

                    # Move to next sibling
                    current = current.next_sibling

                    # Stop if we hit another header
                    if hasattr(current, "name") and current.name in ["h1", "h2", "h3"]:
                        break

            # If no specific pattern found, look for em tags generally
            if not genes:
                em_tags = soup.find_all("em")
                for em in em_tags:
                    text = em.get_text(strip=True)
                    if text:
                        # Check if this is a comma-separated list
                        if "," in text and len(text) > 20:
                            # This might be a gene list
                            potential_genes = [g.strip() for g in text.split(",")]
                            # Only process if it looks like a gene list (multiple valid genes)
                            valid_genes = []
                            for gene in potential_genes:
                                cleaned_gene = self.clean_gene_symbol(gene)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene,
                                ):
                                    valid_genes.append(cleaned_gene)
                            # If we found multiple valid genes, add them all
                            if (
                                len(valid_genes) > 5
                            ):  # Threshold to confirm it's a gene list
                                genes.extend(valid_genes)
                        elif len(text) <= 20:
                            # Single gene
                            cleaned_gene = self.clean_gene_symbol(text)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # Look for gene information in other elements
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li", "strong"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower()
                        for keyword in ["gene", "panel", "exome"]
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

            logger.info(f"Extracted {len(unique_genes)} genes from CEGAT panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing CEGAT panel: {e}")
            raise
