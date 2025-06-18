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
            # Make request with specific network error handling
            try:
                response = requests.get(self.url, timeout=30)
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(
                    f"Network request to Blueprint URL failed: {self.url} - {e}"
                )
                raise  # Re-raise the exception to be caught by the master runner

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")

            genes = []

            # Look for specific gene panel content areas
            panel_content_found = False

            # Strategy 1: Look for "Panel Content" or "Genes" sections
            for keyword in [
                "panel content",
                "genes included",
                "gene list",
                "genes tested",
            ]:
                elements = soup.find_all(
                    string=lambda text, kw=keyword: text and kw in text.lower()
                )
                for element in elements:
                    parent = element.parent
                    logger.info(f"Found '{keyword}' text in element: {parent.name}")

                    # Look for table or list after this element
                    current = parent
                    search_depth = 0
                    while current and search_depth < 10:
                        # Check for tables
                        next_table = (
                            current.find_next("table")
                            if hasattr(current, "find_next")
                            else None
                        )
                        if next_table:
                            logger.info(f"Found table after '{keyword}' text")
                            genes.extend(self._extract_genes_from_table(next_table))
                            panel_content_found = True
                            break

                        # Check for lists
                        next_list = (
                            current.find_next(["ul", "ol"])
                            if hasattr(current, "find_next")
                            else None
                        )
                        if next_list:
                            logger.info(f"Found list after '{keyword}' text")
                            for li in next_list.find_all("li"):
                                gene_text = li.get_text(strip=True)
                                genes.extend(self._extract_genes_from_text(gene_text))
                            panel_content_found = True
                            break

                        current = getattr(current, "next_sibling", None)
                        search_depth += 1

                    if panel_content_found:
                        break

                if panel_content_found:
                    break

            # Strategy 2: Look for tables with gene-like content
            if not genes:
                for table in soup.find_all("table"):
                    # Check if table likely contains genes
                    table_text = table.get_text().upper()
                    if any(
                        indicator in table_text
                        for indicator in ["GENE", "SYMBOL", "OMIM", "TRANSCRIPT"]
                    ):
                        logger.info("Found table with gene indicators")
                        genes.extend(self._extract_genes_from_table(table))

            # Strategy 3: Look for specific div containers with gene lists
            if not genes:
                gene_containers = soup.find_all(
                    "div",
                    class_=lambda x: x
                    and any(
                        term in str(x).lower() for term in ["gene", "panel", "content"]
                    ),
                )
                for container in gene_containers:
                    container_genes = []
                    for element in container.find_all(["span", "p", "li", "td"]):
                        text = element.get_text(strip=True)
                        container_genes.extend(self._extract_genes_from_text(text))
                    if len(container_genes) > 10:  # Likely found a gene list
                        genes.extend(container_genes)
                        break

            # Remove duplicates while preserving order
            unique_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)

            # Final validation - ensure we got real genes
            if len(unique_genes) > 0:
                logger.info(
                    f"Extracted {len(unique_genes)} genes from Blueprint Genetics panel"
                )
            else:
                logger.warning("No genes found, attempting fallback broad scan")
                # Fallback to broader scan if needed
                for element in soup.find_all(["td", "li", "span"]):
                    text = element.get_text(strip=True)
                    if 2 <= len(text) <= 10 and text.isupper():
                        cleaned = self.clean_gene_symbol(text)
                        if (
                            cleaned
                            and self.validate_gene_symbol(cleaned)
                            and self._is_likely_gene(cleaned)
                        ):
                            unique_genes.append(cleaned)
                            seen.add(cleaned)

            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Blueprint Genetics panel: {e}")
            raise

    def _extract_genes_from_table(self, table) -> list[str]:
        """Extract genes from a table element."""
        genes = []
        cells = table.find_all(["td", "th"])
        for cell in cells:
            text = cell.get_text(strip=True)
            genes.extend(self._extract_genes_from_text(text))
        return genes

    def _extract_genes_from_text(self, text: str) -> list[str]:
        """Extract potential genes from a text string."""
        genes = []
        if not text or len(text) > 100:  # Skip very long text
            return genes

        # Handle comma-separated lists
        if "," in text:
            parts = text.split(",")
        else:
            parts = [text]

        for part in parts:
            part = part.strip()
            if 2 <= len(part) <= 15:
                cleaned = self.clean_gene_symbol(part)
                if (
                    cleaned
                    and self.validate_gene_symbol(cleaned)
                    and self._is_likely_gene(cleaned)
                ):
                    genes.append(cleaned)

        return genes

    def _is_likely_gene(self, text: str) -> bool:
        """
        Additional validation to check if text is likely a gene symbol.
        Returns False for known non-gene terms.
        """
        # List of terms that are definitely not genes
        non_gene_terms = {
            # UI/Navigation elements
            "MENU",
            "SEARCH",
            "LOGIN",
            "LOGOUT",
            "SUBMIT",
            "CANCEL",
            "CLOSE",
            "NEXT",
            "PREV",
            "BACK",
            "HOME",
            "ABOUT",
            "CONTACT",
            "HELP",
            "FAQ",
            "NEWS",
            "EVENTS",
            "CAREERS",
            "PRICING",
            "CART",
            "CHECKOUT",
            "PAY",
            # Form elements
            "FIRST",
            "LAST",
            "NAME",
            "EMAIL",
            "PHONE",
            "ADDRESS",
            "CITY",
            "STATE",
            "ZIP",
            "COUNTRY",
            "USERNAME",
            "PASSWORD",
            "CONFIRM",
            "CAPTCHA",
            # Medical terms that aren't genes
            "TEST",
            "TESTS",
            "TESTING",
            "PANEL",
            "PANELS",
            "ANALYSIS",
            "REPORT",
            "RESULTS",
            "CLINICAL",
            "MEDICAL",
            "DIAGNOSTIC",
            "SCREENING",
            # Database/Reference terms
            "OMIM",
            "HGMD",
            "CLINVAR",
            "REFSEQ",
            "GENBANK",
            "PUBMED",
            "DOI",
            # Company/Legal terms
            "LLC",
            "INC",
            "LTD",
            "CORP",
            "COPYRIGHT",
            "RIGHTS",
            "RESERVED",
            "PRIVACY",
            "TERMS",
            "CONDITIONS",
            "POLICY",
            "DISCLAIMER",
            # Technical terms
            "VERSION",
            "UPDATE",
            "DOWNLOAD",
            "UPLOAD",
            "PDF",
            "CSV",
            "JSON",
            "API",
            "SDK",
            "URL",
            "URI",
            "HTTP",
            "HTTPS",
            "SSL",
            "TLS",
        }

        # Check if it's a known non-gene term
        if text.upper() in non_gene_terms:
            return False

        # Check for patterns that indicate non-genes
        # RS numbers (dbSNP identifiers)
        if text.startswith("RS") and len(text) > 3 and text[2:].isdigit():
            return False

        # Concatenated words (like FIRSTLAST, MENUPATIENTS)
        suspicious_patterns = [
            "MENU",
            "FIRST",
            "LAST",
            "ADD",
            "REMOVE",
            "SEARCH",
            "NEWS",
            "CAREER",
            "CONTACT",
            "PRICING",
            "PATIENT",
        ]
        for pattern in suspicious_patterns:
            if pattern in text and len(text) > len(pattern) + 2:
                return False

        # Medical specialties (when they appear as standalone terms)
        specialties = [
            "CARDIOLOGY",
            "DERMATOLOGY",
            "ENDOCRINOLOGY",
            "GASTROENTEROLOGY",
            "HEMATOLOGY",
            "IMMUNOLOGY",
            "NEPHROLOGY",
            "NEUROLOGY",
            "ONCOLOGY",
            "OPHTHALMOLOGY",
            "PEDIATRICS",
            "PSYCHIATRY",
            "PULMONOLOGY",
            "RADIOLOGY",
        ]
        if text in specialties:
            return False

        # Disease/condition terms that are not genes
        diseases = [
            "ANEMIA",
            "CARCINOMA",
            "MELANOMA",
            "LYMPHOMA",
            "NEUTROPENIA",
            "THROMBOCYTOSIS",
            "SPHEROCYTOSIS",
            "OSTEOPETROSIS",
            "HEMOCHROMATOSIS",
            "ALBINISM",
            "ANGIOEDEMA",
            "CYLINDROMATOSIS",
            "DEAFNESS",
            "ELLIPSOCYTOSIS",
            "ERYTHROCYTOSIS",
            "ERTHYROCYTOSIS",
            "ERYTHROPOIETIC",
            "FIBROMATOSIS",
            "MYOPATHY",
            "MYXOMA",
            "NYSTAGMUS",
            "OLIGONDONTIA",
            "PERIODONTITIS",
            "PIEBALDISM",
            "PNEUMOTHORAX",
            "POLYPOSIS",
            "PROTOPORPHYRIA",
            "SCHWANNOMATOSIS",
            "SCOLIOSIS",
            "SEIZURES",
            "SPASTICITY",
            "WARTS",
            "SITOSTEROLEMIA",
            "THROMBOPHILIA",
            "THROMBOCYTHEMIA",
            "NEUTROPHILIA",
            "ACRODYSOSTOSIS",
            "AFIBRINOGENEMIA",
            "HEMOGLOBIN",
        ]
        if text in diseases:
            return False

        # Descriptive terms
        descriptive = [
            "FAMILIAL",
            "JUVENILE",
            "ADULT-ONSET",
            "CHILDHOOD-ONSET",
            "CONGENITAL",
            "MULTIPLE",
            "ISOLATED",
            "ATYPICAL",
            "INTRACARDIAC",
            "INHERITANCE",
            "FINNISH",
            "NORWEGIAN",
            "PULMONARY",
            "PAPILLARY",
            "PARATHYROID",
            "PHOTOSENSITIVE",
            "PLATELET-TYPE",
            "RECURRENT",
            "INFECTIONS",
            "MULTISYSTEM",
            "NON-HODGKIN",
            "GINGIVAL",
            "OCULOCUTANEOUS",
            "PARAGANGLIOMAS",
            "LYMPEHEDEMA",
            "VKORC1-RELATED",
        ]
        if text in descriptive:
            return False

        # Special cases - these might look like genes but aren't
        if text in ["AD", "AR", "XL"]:  # These could be inheritance patterns, not genes
            return False

        # If it made it this far, it's likely a gene
        return True
