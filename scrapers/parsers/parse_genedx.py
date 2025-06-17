"""
GeneDx parser.

This module extracts gene information from GeneDx panel pages.
Uses Selenium since the content is loaded via JavaScript.
"""

import logging
import os
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
            # Temporarily set no_proxy environment variable for localhost
            original_no_proxy = os.environ.get("no_proxy", "")
            os.environ["no_proxy"] = "localhost,127.0.0.1"

            # Setup Chrome options for headless browsing (WSL-compatible)
            chrome_options = Options()
            chrome_options.add_argument("--headless")
            chrome_options.add_argument("--no-sandbox")
            chrome_options.add_argument("--disable-dev-shm-usage")
            chrome_options.add_argument("--disable-gpu")
            chrome_options.add_argument("--disable-web-security")
            chrome_options.add_argument("--disable-features=VizDisplayCompositor")
            chrome_options.add_argument("--window-size=1920,1080")
            chrome_options.add_argument("--remote-debugging-port=9223")
            chrome_options.add_argument("--disable-background-timer-throttling")
            chrome_options.add_argument("--disable-renderer-backgrounding")
            chrome_options.add_argument("--disable-backgrounding-occluded-windows")
            # Configure proxy for Charite network
            chrome_options.add_argument("--proxy-server=http://proxy.charite.de:8080")
            chrome_options.add_argument("--proxy-bypass-list=localhost,127.0.0.1")
            chrome_options.add_argument(
                "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
            )

            # Initialize driver - try manual ChromeDriver first, then fallback
            try:
                # Try manually downloaded ChromeDriver for Chrome 137
                service = Service("/tmp/chromedriver-linux64/chromedriver")
                driver = webdriver.Chrome(service=service, options=chrome_options)
            except Exception:
                try:
                    # Fallback to WebDriver Manager
                    service = Service(ChromeDriverManager().install())
                    driver = webdriver.Chrome(service=service, options=chrome_options)
                except Exception as e:
                    logger.error(f"Failed to initialize ChromeDriver: {e}")
                    raise

            # Load the page
            driver.get(self.url)

            # Wait for content to load
            WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.TAG_NAME, "body"))
            )

            # Additional wait for dynamic content
            time.sleep(3)

            genes = []

            # Wait for the custom panel interface to load
            try:
                WebDriverWait(driver, 15).until(
                    EC.presence_of_element_located(
                        (By.CSS_SELECTOR, ".list-item__title")
                    )
                )
                time.sleep(2)  # Additional wait for dynamic content
                logger.info("Custom panel interface loaded")
            except Exception:
                logger.warning(
                    "Custom panel interface not found, trying alternative selectors"
                )

            # Extract genes from all pages of the custom panel checkbox list
            page_num = 1
            max_pages = 10  # Safety limit to prevent infinite loops

            while page_num <= max_pages:
                logger.info(f"Processing page {page_num}")

                # Extract genes from current page
                gene_title_elements = driver.find_elements(
                    By.CSS_SELECTOR, ".list-item__title"
                )
                current_page_genes = 0

                for element in gene_title_elements:
                    text = element.text.strip()
                    if text:
                        # R script removes " (GREM1)" pattern and other annotations
                        text = text.replace(" (GREM1)", "").replace(
                            "SCG5 (GREM1)", "GREM1"
                        )
                        cleaned_gene = self.clean_gene_symbol(text)
                        if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                            if cleaned_gene not in genes:  # Avoid duplicates
                                genes.append(cleaned_gene)
                                current_page_genes += 1

                logger.info(f"Found {current_page_genes} new genes on page {page_num}")

                # Try to find and click the next page button
                try:
                    # Look for next page button that's not disabled
                    next_buttons = driver.find_elements(
                        By.CSS_SELECTOR, "button[aria-label='']:not([disabled])"
                    )
                    next_button = None

                    for btn in next_buttons:
                        # Check if it contains a right arrow icon
                        arrow_icons = btn.find_elements(
                            By.CSS_SELECTOR, "svg path[d*='M7.293 14.707']"
                        )
                        if arrow_icons:
                            next_button = btn
                            break

                    if next_button:
                        logger.info("Clicking next page button")
                        driver.execute_script("arguments[0].click();", next_button)
                        time.sleep(3)  # Wait for page to load
                        page_num += 1
                    else:
                        logger.info("No more pages available")
                        break

                except Exception as e:
                    logger.info(f"No more pages or error navigating: {e}")
                    break

            # If no gene titles found, try checkbox labels
            if not genes:
                checkbox_labels = driver.find_elements(
                    By.CSS_SELECTOR, "input[id*='selection-'] + span"
                )
                for element in checkbox_labels:
                    text = element.text.strip()
                    if text:
                        cleaned_gene = self.clean_gene_symbol(text)
                        if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                            genes.append(cleaned_gene)

            # Try parsing from JavaScript data (genes attribute)
            if not genes:
                try:
                    # Look for genes in data attributes
                    custom_panel_element = driver.find_element(
                        By.CSS_SELECTOR, "[genes]"
                    )
                    genes_attribute = custom_panel_element.get_attribute("genes")
                    if genes_attribute:
                        import json

                        gene_list = json.loads(genes_attribute)
                        for gene in gene_list:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                except Exception:
                    pass

            # If no genes found, try other selectors
            if not genes:
                selectors = [
                    "span.gene-name",  # Original R script selector
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
            # Restore original no_proxy environment variable
            if "original_no_proxy" in locals():
                if original_no_proxy:
                    os.environ["no_proxy"] = original_no_proxy
                else:
                    os.environ.pop("no_proxy", None)
            if driver:
                driver.quit()
