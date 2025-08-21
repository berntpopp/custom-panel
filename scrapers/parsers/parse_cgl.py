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

            # Look for any text containing patterns like GENE_NM_XXXXXX
            # We'll search the entire page content for these patterns
            page_text = soup.get_text()

            import re

            # Find all patterns that look like GENE_NM_XXXXXX or GENE_TRANSCRIPT
            # Changed {2,10} to {1,10} to allow 2-character genes like FH
            gene_patterns = re.findall(r"\b([A-Z][A-Z0-9]{1,10})_[A-Z0-9_]+", page_text)

            if gene_patterns:
                logger.info(f"Found {len(gene_patterns)} gene patterns in page")
                for gene_symbol in gene_patterns:
                    cleaned_gene = self.clean_gene_symbol(gene_symbol)
                    if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                        genes.append(cleaned_gene)

            # If no underscore patterns found, look for GENES TARGETED sections specifically
            if not genes:
                # Look for elements with "genes targeted" and extract following content
                for element in soup.find_all(
                    string=lambda text: text and "genes targeted" in text.lower(),
                ):
                    parent = element.parent
                    logger.info(
                        "Found 'GENES TARGETED' text, looking for gene content...",
                    )

                    # Look for content in following elements
                    for next_element in parent.find_all_next(
                        ["p", "div", "li", "span"], limit=20,
                    ):
                        text = next_element.get_text(strip=True)
                        if not text:
                            continue

                        # Skip navigation and descriptive text
                        if any(
                            skip in text.lower()
                            for skip in [
                                "table of contents",
                                "toggle",
                                "overview",
                                "test requirements",
                                "turn-around time",
                                "results reporting",
                                "method",
                                "built with",
                                "copyright",
                                "cancer genetics and genomics laboratory",
                            ]
                        ):
                            continue

                        # Look for semicolon-separated gene lists or gene patterns
                        if (
                            ";" in text and len(text) < 1000
                        ):  # Reasonable length for gene list
                            parts = text.split(";")
                            for part in parts:
                                part = part.strip()
                                if "_" in part:
                                    gene_symbol = part.split("_")[0].strip()
                                    if gene_symbol and len(gene_symbol) >= 3:
                                        cleaned_gene = self.clean_gene_symbol(
                                            gene_symbol,
                                        )
                                        if cleaned_gene and self.validate_gene_symbol(
                                            cleaned_gene,
                                        ):
                                            genes.append(cleaned_gene)

                        # Stop after finding some genes
                        if len(genes) > 5:
                            break

                    if genes:
                        break

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
