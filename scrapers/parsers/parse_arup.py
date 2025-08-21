"""
ARUP Laboratories parser.

This module extracts gene information from ARUP panel pages.
Looks for genes in tables following elements with id="genes-tested1".
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class ArupParser(BaseParser):
    """Parser for ARUP Laboratories panels."""

    def parse(self) -> list[str]:
        """
        Parse ARUP panel page and extract gene symbols.

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
                logger.error(f"Network request to ARUP URL failed: {self.url} - {e}")
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Look for "Genes Tested" text and find the following table
            genes_tested_found = False

            # Try multiple approaches to find "Genes Tested"
            for element in soup.find_all(
                string=lambda text: text and "genes tested" in text.lower(),
            ):
                parent = element.parent
                logger.info(f"Found 'Genes Tested' text in element: {parent.name}")

                # Look for ALL tables after this element, not just the first one
                current = parent
                tables_found = 0

                # Find all tables that come after "Genes Tested"
                while current and tables_found < 5:  # Limit to avoid infinite loops
                    next_table = (
                        current.find_next("table")
                        if hasattr(current, "find_next")
                        else None
                    )
                    if next_table:
                        tables_found += 1
                        logger.info(
                            f"Found table {tables_found} after 'Genes Tested' text",
                        )

                        # Extract genes from table cells
                        cells = next_table.find_all(["td", "th"])
                        table_genes = []
                        for cell in cells:
                            text = cell.get_text(strip=True)
                            # Gene symbols are typically 3-10 characters
                            if text and 3 <= len(text) <= 10:
                                # Skip common false positives
                                false_positives = {
                                    "GENE",
                                    "GENES",
                                    "TEST",
                                    "PANEL",
                                    "CANCER",
                                    "ADB",
                                    "ADC",
                                    "BREAST",
                                    "BREASTA",
                                    "OVARIAN",
                                    "COLON",
                                    "LUNG",
                                    "PROSTATE",
                                    "KIDNEY",
                                    "LIVER",
                                    "BRAIN",
                                    "TUMOR",
                                    "TUMORS",
                                    "SNVS",
                                    "GIST",
                                }
                                if text.upper() in false_positives:
                                    continue

                                # Check if it looks like a gene symbol (mostly letters/numbers)
                                if (
                                    text.replace("1", "")
                                    .replace("2", "")
                                    .replace("3", "")
                                    .replace("4", "")
                                    .replace("5", "")
                                    .isalpha()
                                ):
                                    cleaned_gene = self.clean_gene_symbol(text)
                                    if cleaned_gene and self.validate_gene_symbol(
                                        cleaned_gene,
                                    ):
                                        table_genes.append(cleaned_gene)

                        logger.info(
                            f"Table {tables_found} contributed {len(table_genes)} genes",
                        )
                        genes.extend(table_genes)

                        # Move current to the table we just processed to find the next one
                        current = next_table
                    else:
                        break

                if genes:
                    genes_tested_found = True
                    break

                if genes_tested_found:
                    break

            # Alternative approach: look for element with id containing "genes-tested"
            if not genes:
                for element in soup.find_all(
                    attrs={"id": lambda x: x and "genes-tested" in x.lower()},
                ):
                    logger.info(
                        f"Found element with genes-tested id: {element.get('id')}",
                    )

                    # Look for table after this element
                    table = element.find_next("table")
                    if table:
                        logger.info("Found table after genes-tested element")
                        cells = table.find_all(["td", "th"])
                        for cell in cells:
                            text = cell.get_text(strip=True)
                            if text and 3 <= len(text) <= 10:
                                # Skip common false positives
                                false_positives = {
                                    "GENE",
                                    "GENES",
                                    "TEST",
                                    "PANEL",
                                    "CANCER",
                                    "ADB",
                                    "ADC",
                                    "BREAST",
                                    "BREASTA",
                                    "OVARIAN",
                                    "COLON",
                                    "LUNG",
                                    "PROSTATE",
                                    "KIDNEY",
                                    "LIVER",
                                    "BRAIN",
                                    "TUMOR",
                                    "TUMORS",
                                    "SNVS",
                                    "GIST",
                                }
                                if text.upper() in false_positives:
                                    continue

                                # Check if it looks like a gene symbol (mostly letters/numbers)
                                if (
                                    text.replace("1", "")
                                    .replace("2", "")
                                    .replace("3", "")
                                    .replace("4", "")
                                    .replace("5", "")
                                    .isalpha()
                                ):
                                    cleaned_gene = self.clean_gene_symbol(text)
                                    if cleaned_gene and self.validate_gene_symbol(
                                        cleaned_gene,
                                    ):
                                        genes.append(cleaned_gene)
                        break

            # Remove duplicates while preserving order
            unique_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)

            logger.info(f"Extracted {len(unique_genes)} genes from ARUP panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing ARUP panel: {e}")
            raise
