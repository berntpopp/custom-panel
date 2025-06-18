"""
Abstract base class for all commercial panel parsers.

This module defines the BaseParser class that all individual parsers must inherit from.
"""

from abc import ABC, abstractmethod
from typing import Any


class BaseParser(ABC):
    """
    Abstract base class for commercial diagnostic panel parsers.

    All parser implementations must inherit from this class and implement the parse method.
    """

    def __init__(self, url: str, config: dict[str, Any] | None = None) -> None:
        """
        Initialize the parser with URL and optional configuration.

        Args:
            url: The URL to scrape data from
            config: Optional configuration dictionary
        """
        self.url = url
        self.config = config or {}

    @abstractmethod
    def parse(self) -> list[str]:
        """
        Parse the target URL and extract gene symbols.

        Returns:
            List of gene symbols found on the page

        Raises:
            Exception: If parsing fails
        """
        pass

    def clean_gene_symbol(self, gene: str) -> str:
        """
        Clean and standardize a gene symbol.

        Args:
            gene: Raw gene symbol string

        Returns:
            Cleaned gene symbol
        """
        if not gene:
            return ""

        # Remove common artifacts
        gene = gene.strip()
        gene = gene.replace("*", "")  # Remove asterisks
        gene = gene.replace("#", "")  # Remove hash symbols
        gene = gene.split("(")[0].strip()  # Remove parenthetical content
        gene = gene.split(",")[0].strip()  # Take first gene if comma-separated

        # Convert to uppercase for consistency
        return gene.upper()

    def validate_gene_symbol(self, gene: str) -> bool:
        """
        Validate that a string looks like a valid gene symbol.

        Args:
            gene: Gene symbol to validate

        Returns:
            True if the gene symbol appears valid
        """
        if not gene:
            return False

        # Basic validation: should be alphanumeric with possible hyphens/underscores
        # and between 1-20 characters
        if len(gene) < 1 or len(gene) > 20:
            return False

        # Should contain at least one letter
        if not any(c.isalpha() for c in gene):
            return False

        # Should only contain valid characters
        valid_chars = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_")
        if not all(c in valid_chars for c in gene):
            return False

        # Skip common non-gene terms
        skip_terms = {
            "GENE",
            "DNA",
            "RNA",
            "PANEL",
            "TEST",
            "CANCER",
            "TUMOR",
            "DISEASE",
            "SYNDROME",
            "MUTATION",
            "VARIANT",
            "SEQUENCE",
            "INCLUDES",
            "GENES",
            "FOR",
            "HEREDITARY",
            "RISK",
            "ASSESSMENT",
            "COMPREHENSIVE",
            "FULL",
        }
        if gene in skip_terms:
            return False

        return True
