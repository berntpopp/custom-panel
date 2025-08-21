"""
Centogene parser.

This module extracts gene information from Centogene panel pages.
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


class CentogeneParser(BaseParser):
    """Parser for Centogene panels (requires Selenium for JavaScript content)."""

    def parse(self) -> list[str]:
        """
        Parse Centogene panel page and extract gene symbols.

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
            chrome_options.add_argument("--remote-debugging-port=9224")
            chrome_options.add_argument("--disable-background-timer-throttling")
            chrome_options.add_argument("--disable-renderer-backgrounding")
            chrome_options.add_argument("--disable-backgrounding-occluded-windows")
            # Configure proxy for Charite network
            chrome_options.add_argument("--proxy-server=http://proxy.charite.de:8080")
            chrome_options.add_argument("--proxy-bypass-list=localhost,127.0.0.1")
            chrome_options.add_argument(
                "--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
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
                EC.presence_of_element_located((By.TAG_NAME, "body")),
            )

            # Additional wait for dynamic content
            time.sleep(3)

            genes = []

            # Look for the "Show all genes" expansion panel
            try:
                WebDriverWait(driver, 15).until(
                    EC.presence_of_element_located(
                        (By.CSS_SELECTOR, ".mat-expansion-panel"),
                    ),
                )
                logger.info("Found expansion panel")

                # Find expansion panel with "Show all genes" text
                expansion_panels = driver.find_elements(
                    By.CSS_SELECTOR,
                    ".mat-expansion-panel",
                )
                genes_panel = None

                for panel in expansion_panels:
                    panel_text = panel.text
                    if "show all genes" in panel_text.lower():
                        genes_panel = panel
                        break

                if genes_panel:
                    logger.info("Found 'Show all genes' panel")

                    # Try to click the expansion panel header to expand it
                    try:
                        header = genes_panel.find_element(
                            By.CSS_SELECTOR,
                            ".mat-expansion-panel-header",
                        )
                        driver.execute_script("arguments[0].click();", header)
                        time.sleep(2)
                        logger.info("Clicked to expand panel")
                    except Exception as e:
                        logger.warning(f"Could not click expansion panel: {e}")

                    # Look for gene content in the expanded panel
                    panel_body = genes_panel.find_element(
                        By.CSS_SELECTOR,
                        ".mat-expansion-panel-body",
                    )
                    if panel_body:
                        gene_text = panel_body.text.strip()
                        logger.info(f"Found gene text: {gene_text[:100]}...")

                        # Parse comma-separated genes
                        if gene_text:
                            gene_parts = gene_text.split(",")
                            for gene_part in gene_parts:
                                gene_part = gene_part.strip()
                                cleaned_gene = self.clean_gene_symbol(gene_part)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene,
                                ):
                                    genes.append(cleaned_gene)

            except Exception as e:
                logger.warning(f"Could not find expansion panel: {e}")

            # Fallback: look for genes in any visible text on the page
            if not genes:
                logger.info("Trying fallback approach to find genes")
                page_text = driver.find_element(By.TAG_NAME, "body").text
                import re

                # Look for gene-like patterns in the entire page
                potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,7}\b", page_text)

                # Filter potential genes to only include likely gene symbols
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

            logger.info(f"Extracted {len(unique_genes)} genes from Centogene panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Centogene panel: {e}")
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
