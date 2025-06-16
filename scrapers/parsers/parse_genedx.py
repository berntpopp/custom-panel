"""
GeneDx parser.

This module extracts gene information from GeneDx panel pages.
Uses Selenium since the content is loaded via JavaScript.
"""

import logging
import time

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class GenedxParser(BaseParser):
    """Parser for GeneDx panels (requires Selenium for JavaScript content)."""

    def parse(self) -> list[str]:
        """
        Parse GeneDx panel page and extract gene symbols.

        Returns:
            List of gene symbols found on the page

        Raises:
            Exception: If parsing fails
        """
        driver = None
        try:
            # Setup Chrome options for headless browsing
            chrome_options = Options()
            chrome_options.add_argument("--headless")
            chrome_options.add_argument("--no-sandbox")
            chrome_options.add_argument("--disable-dev-shm-usage")
            chrome_options.add_argument("--disable-gpu")
            chrome_options.add_argument("--window-size=1920,1080")
            chrome_options.add_argument(
                "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
            )

            # Initialize driver
            service = Service(ChromeDriverManager().install())
            driver = webdriver.Chrome(service=service, options=chrome_options)

            # Load the page
            driver.get(self.url)

            # Wait for content to load
            WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.TAG_NAME, "body"))
            )

            # Additional wait for dynamic content
            time.sleep(3)

            genes = []

            # Find spans with class="gene-name"
            gene_name_elements = driver.find_elements(By.CSS_SELECTOR, "span.gene-name")
            for element in gene_name_elements:
                text = element.text.strip()
                if text:
                    cleaned_gene = self.clean_gene_symbol(text)
                    if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                        genes.append(cleaned_gene)

            # If no gene-name elements found, try other selectors
            if not genes:
                selectors = [
                    ".gene",
                    "[data-gene]",
                    "span[class*='gene']",
                    "td.gene",
                    ".gene-symbol",
                ]

                for selector in selectors:
                    try:
                        elements = driver.find_elements(By.CSS_SELECTOR, selector)
                        for element in elements:
                            text = element.text.strip()
                            if text and len(text) <= 20:
                                cleaned_gene = self.clean_gene_symbol(text)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene
                                ):
                                    genes.append(cleaned_gene)
                    except Exception:
                        continue

            # Look for genes in tables
            if not genes:
                try:
                    table_cells = driver.find_elements(By.TAG_NAME, "td")
                    for cell in table_cells:
                        text = cell.text.strip()
                        if text and len(text) <= 20:
                            cleaned_gene = self.clean_gene_symbol(text)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                except Exception:
                    pass

            # Look for genes in JavaScript data
            if not genes:
                scripts = driver.find_elements(By.TAG_NAME, "script")
                for script in scripts:
                    script_content = script.get_attribute("innerHTML")
                    if script_content and "gene" in script_content.lower():
                        import re

                        potential_genes = re.findall(
                            r"\b[A-Z][A-Z0-9]{1,7}\b", script_content
                        )
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

            logger.info(f"Extracted {len(unique_genes)} genes from GeneDx panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing GeneDx panel: {e}")
            raise
        finally:
            if driver:
                driver.quit()
