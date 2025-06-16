"""
Prevention Genetics parser.

This module extracts gene information from Prevention Genetics panel pages.
Looks for genes in tables with id="genes-div".
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class PreventionParser(BaseParser):
    """Parser for Prevention Genetics panels."""

    def parse(self) -> list[str]:
        """
        Parse Prevention Genetics panel page and extract gene symbols.

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
                    f"Network request to Prevention Genetics URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Find table with id="genes-div"
            genes_table = soup.find(id="genes-div")
            if genes_table:
                # Look for 'Official Gene Symbol' column
                table = (
                    genes_table.find("table")
                    if genes_table.name != "table"
                    else genes_table
                )
                if table:
                    headers = table.find_all("th")
                    gene_col_idx = None

                    for i, header in enumerate(headers):
                        if "official gene symbol" in header.get_text().lower():
                            gene_col_idx = i
                            break

                    if gene_col_idx is not None:
                        rows = table.find_all("tr")
                        for row in rows[1:]:  # Skip header row
                            cells = row.find_all(["td", "th"])
                            if len(cells) > gene_col_idx:
                                gene_text = cells[gene_col_idx].get_text(strip=True)
                                if gene_text:
                                    cleaned_gene = self.clean_gene_symbol(gene_text)
                                    if cleaned_gene and self.validate_gene_symbol(
                                        cleaned_gene
                                    ):
                                        genes.append(cleaned_gene)

            # If specific table not found, look for genes generally
            if not genes:
                tables = soup.find_all("table")
                for table in tables:
                    cells = table.find_all(["td", "th"])
                    for cell in cells:
                        text = cell.get_text(strip=True)
                        if text and len(text) <= 20:
                            cleaned_gene = self.clean_gene_symbol(text)
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
                f"Extracted {len(unique_genes)} genes from Prevention Genetics panel"
            )
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Prevention Genetics panel: {e}")
            raise
