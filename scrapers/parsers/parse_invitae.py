"""
Invitae parser.

This module extracts gene information from Invitae panel pages.
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


class InvitaeParser(BaseParser):
    """Parser for Invitae panels (requires Selenium for JavaScript content)."""

    def parse(self) -> list[str]:
        """
        Parse Invitae panel page and extract gene symbols.

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

            # R script uses: //meta[contains(@content,"AIP")]
            # Look for meta tags containing "AIP" in their content
            meta_elements = driver.find_elements(By.TAG_NAME, "meta")
            for meta in meta_elements:
                content = meta.get_attribute("content") or ""
                if "AIP" in content:
                    import re

                    # Split content and extract gene symbols
                    potential_genes = re.split(r"[,;\s]+", content)
                    for gene in potential_genes:
                        gene = gene.strip()
                        if gene:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                    break

            # Look for genes in JavaScript variables
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

            # Look for genes in page content
            if not genes:
                # Try to find gene lists in various elements
                selectors = [
                    "span.gene-name",
                    ".gene",
                    "[data-gene]",
                    "td",
                    "li",
                    "div",
                    "p",
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

            # Remove duplicates while preserving order
            unique_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)

            logger.info(f"Extracted {len(unique_genes)} genes from Invitae panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Invitae panel: {e}")
            raise
        finally:
            if driver:
                driver.quit()
