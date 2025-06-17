"""
CGL (Cancer Genetics Lab) parser.

This module extracts gene information from CGL panel pages.
Looks for genes in paragraphs following "GENES TARGETED" headers.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class CglParser(BaseParser):
    """Parser for CGL (Cancer Genetics Lab) panels."""

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
            try:
                response = requests.get(self.url, timeout=30)
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(f"Network request to CGL URL failed: {self.url} - {e}")
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # R script uses: //h2[contains(text(),"GENES TARGETED")]//following::p
            # Find h2 with "GENES TARGETED" text and look for following p tags
            genes_targeted_headers = soup.find_all(
                "h2", string=lambda text: text and "genes targeted" in text.lower()
            )

            for header in genes_targeted_headers:
                # Find all following paragraphs, not just the next one
                for p_tag in header.find_all_next("p"):
                    text = p_tag.get_text(strip=True)
                    if text:
                        # R script applies several filters and cleaning
                        if any(
                            skip in text
                            for skip in [
                                "Single nucleotide variants",
                                "Variants not presumed by nature",
                                "Genomic DNA is subjected to FFPE repair",
                                "Variants are validated",
                            ]
                        ):
                            continue

                        # Remove R script patterns
                        text = text.replace("Entire coding region: ", "")
                        text = text.replace("Partial Genes: ", "")
                        text = text.replace(
                            "Copy Number Variants (not applicable to FFPE specimens): Copy number variants are called in all of the above genes as well as:",
                            "",
                        )

                        # Split by semicolon like R script
                        parts = text.split(";")
                        for part in parts:
                            part = part.strip()
                            if part:
                                # Remove underscore patterns like R script
                                if "_" in part:
                                    part = part.split("_")[0]
                                cleaned_gene = self.clean_gene_symbol(part)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene
                                ):
                                    genes.append(cleaned_gene)

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
