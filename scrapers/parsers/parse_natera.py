"""
Natera parser.

This module extracts gene information from Natera panel pages.
Uses Selenium as the content is dynamic and looks for text following "Genes:" labels.
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


class NateraParser(BaseParser):
    """Parser for Natera panels (requires Selenium for JavaScript content)."""

    def parse(self) -> list[str]:
        """
        Parse Natera panel page and extract gene symbols.

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

            # Find strong tags with "Genes:" and extract following text
            strong_elements = driver.find_elements(By.TAG_NAME, "strong")
            for strong in strong_elements:
                if "genes:" in strong.text.lower():
                    # Get the parent element and extract text after "Genes:"
                    parent = strong.find_element(By.XPATH, "..")
                    text = parent.text

                    # Find the part after "Genes:"
                    if "genes:" in text.lower():
                        gene_part = text.lower().split("genes:")[1].strip()

                        # Parse gene list
                        for gene in gene_part.split(","):
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # If no specific pattern found, look for gene lists generally
            if not genes:
                elements = driver.find_elements(By.CSS_SELECTOR, "div, p, span, li")
                for element in elements:
                    text = element.text
                    if any(
                        keyword in text.lower()
                        for keyword in ["gene", "panel", "ctdna"]
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

            logger.info(f"Extracted {len(unique_genes)} genes from Natera panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Natera panel: {e}")
            raise
        finally:
            if driver:
                driver.quit()
