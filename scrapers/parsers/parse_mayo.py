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
            # Make request with proxy support
            try:
                # Setup proxy configuration for Charite network
                proxies = {
                    "http": "http://proxy.charite.de:8080",
                    "https": "http://proxy.charite.de:8080",
                }
                headers = {
                    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
                }
                response = requests.get(
                    self.url, proxies=proxies, headers=headers, timeout=30
                )
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(
                    f"Network request to Mayo Clinic URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Look for genes in the Genetics Test Information section
            genetics_info = soup.find("div", class_="GeneticsTestInformation-value")
            if genetics_info:
                # Find the paragraph containing the gene list
                paragraphs = genetics_info.find_all("p")
                for paragraph in paragraphs:
                    text = paragraph.get_text(strip=True)
                    # Check if this paragraph contains gene symbols (comma-separated list)
                    if "," in text and any(
                        word.isupper() and len(word) <= 7 for word in text.split()
                    ):
                        logger.info(
                            "Found gene list in Genetics Test Information section"
                        )
                        # Split by comma and extract genes
                        gene_parts = text.split(",")
                        for gene_part in gene_parts:
                            # Clean up annotations like "(including promoters 1A and 1B)"
                            gene_part = gene_part.strip()
                            # Remove parenthetical information
                            import re

                            gene_part = re.sub(r"\s*\([^)]*\)", "", gene_part)

                            cleaned_gene = self.clean_gene_symbol(gene_part)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                        break  # Found the gene list, no need to continue

            # Fallback: Find italic tags containing gene lists (original logic)
            if not genes:
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

            # General search for gene patterns in divs/spans that might contain gene lists
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li"]):
                    text = element.get_text()
                    # Look for comma-separated lists that might be genes
                    if "," in text and len(text) < 2000:  # Reasonable length limit
                        import re

                        # Extract potential gene symbols (2-7 uppercase letters/numbers)
                        potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,6}\b", text)
                        # Only proceed if we find multiple potential genes
                        if len(potential_genes) >= 5:
                            for gene in potential_genes:
                                cleaned_gene = self.clean_gene_symbol(gene)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene
                                ):
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
