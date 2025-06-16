"""
Blueprint Genetics panel parser.

This module extracts gene information from Blueprint Genetics panel pages.
Handles multiple sub-panel URLs and aggregates unique genes.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class BlueprintParser(BaseParser):
    """Parser for Blueprint Genetics panels."""

    def parse(self) -> list[str]:
        """
        Parse Blueprint Genetics panel page and extract gene symbols.

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

            # Find tables with class "table mb-5"
            gene_tables = soup.find_all(
                "table", class_=lambda x: x and "table" in x and "mb-5" in x
            )

            genes = []

            for table in gene_tables:
                # Extract genes from table cells
                cells = table.find_all(["td", "th"])
                for cell in cells:
                    text = cell.get_text(strip=True)
                    if text and len(text) <= 20:  # Likely a gene symbol
                        cleaned_gene = self.clean_gene_symbol(text)
                        if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                            genes.append(cleaned_gene)

            # If no tables found, look for gene lists in other elements
            if not genes:
                # Look for gene listings in various elements
                for element in soup.find_all(["div", "span", "p", "li"]):
                    text = element.get_text(strip=True)

                    # Skip if text is too long (not a gene symbol)
                    if len(text) > 20:
                        continue

                    # Check if it looks like a gene symbol
                    cleaned_gene = self.clean_gene_symbol(text)
                    if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                        genes.append(cleaned_gene)

            # Also look for structured gene data in JSON-LD or data attributes
            scripts = soup.find_all("script", type="application/ld+json")
            for script in scripts:
                try:
                    import json

                    data = json.loads(script.string)
                    # Look for gene mentions in structured data
                    if isinstance(data, dict) and "genes" in str(data).lower():
                        # Extract potential genes from the JSON string
                        import re

                        potential_genes = re.findall(
                            r"\b[A-Z][A-Z0-9]{1,7}\b", str(data)
                        )
                        for gene in potential_genes:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                except Exception:
                    continue

            # Remove duplicates while preserving order
            unique_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)

            logger.info(
                f"Extracted {len(unique_genes)} genes from Blueprint Genetics panel"
            )
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Blueprint Genetics panel: {e}")
            raise
