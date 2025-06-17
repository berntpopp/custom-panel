"""
Natera parser.

This module extracts gene information from Natera panel pages.
Uses Selenium as the content is dynamic and looks for text following "Genes:" labels.
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

            # R script uses: //strong[contains(text(),"Genes:")]/following-sibling::span|//strong[contains(text(),"Genes:")]/following-sibling::text()
            # Find strong tags with "Genes:" and extract following siblings
            try:
                # Try the exact XPath from R script
                gene_elements = driver.find_elements(
                    By.XPATH,
                    '//strong[contains(text(),"Genes:")]/following-sibling::span | //strong[contains(text(),"Genes:")]/following-sibling::text()',
                )
                for element in gene_elements:
                    text = element.text if hasattr(element, "text") else str(element)
                    if text:
                        # Split by comma and parse like R script
                        for gene in text.split(","):
                            gene = gene.strip()
                            if gene:
                                cleaned_gene = self.clean_gene_symbol(gene)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene
                                ):
                                    genes.append(cleaned_gene)
            except Exception:
                # Fallback to original approach
                strong_elements = driver.find_elements(By.TAG_NAME, "strong")
                for strong in strong_elements:
                    if "genes:" in strong.text.lower():
                        # Get following sibling elements
                        try:
                            following_elements = driver.execute_script(
                                """
                                var element = arguments[0];
                                var siblings = [];
                                var sibling = element.nextSibling;
                                while (sibling) {
                                    if (sibling.nodeType === 1 || sibling.nodeType === 3) {
                                        siblings.push(sibling.textContent || sibling.data || '');
                                    }
                                    sibling = sibling.nextSibling;
                                }
                                return siblings.join(' ');
                            """,
                                strong,
                            )

                            if following_elements:
                                for gene in following_elements.split(","):
                                    gene = gene.strip()
                                    if gene:
                                        cleaned_gene = self.clean_gene_symbol(gene)
                                        if cleaned_gene and self.validate_gene_symbol(
                                            cleaned_gene
                                        ):
                                            genes.append(cleaned_gene)
                        except Exception:
                            continue

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
            # Restore original no_proxy environment variable
            if "original_no_proxy" in locals():
                if original_no_proxy:
                    os.environ["no_proxy"] = original_no_proxy
                else:
                    os.environ.pop("no_proxy", None)
            if driver:
                driver.quit()
