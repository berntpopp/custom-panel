"""
Invitae parser.

This module extracts gene information from Invitae panel pages.
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


class InvitaeParser(BaseParser):
    """Parser for Invitae panels (requires Selenium for JavaScript content)."""

    def _is_likely_gene(self, gene: str) -> bool:
        """
        Additional filtering to identify likely gene symbols.

        Args:
            gene: Gene symbol to check

        Returns:
            True if this looks like a real gene symbol
        """
        # Skip common non-gene web terms
        web_terms = {
            "ADMIN",
            "AGENT",
            "AJAX",
            "ANALYTICS",
            "ANGULAR",
            "API",
            "AWS",
            "AZURE",
            "BABEL",
            "BACKBONE",
            "BACKUP",
            "BITBUCKET",
            "BODY",
            "BOOTSTRAP",
            "BUILD",
            "BUNDLE",
            "CACHE",
            "CANSCAN",
            "CAP",
            "CDATA",
            "CD",
            "CDN",
            "CI",
            "CLIA",
            "CMIPECGS",
            "CMMR",
            "CMS",
            "COMPILE",
            "CONFIG",
            "CORS",
            "CRISPR",
            "CRM",
            "CRS",
            "CSS",
            "DB",
            "DEBUG",
            "DELETE",
            "DEPLOY",
            "DEV",
            "DNS",
            "DOCKER",
            "EDTA",
            "EMAIL",
            "EMBER",
            "END",
            "ENV",
            "ERP",
            "ERROR",
            "EVERYONE",
            "FAQ",
            "FETCH",
            "FTP",
            "GCP",
            "GENAF",
            "GET",
            "GITHUB",
            "GIUGTLBJ",
            "GITLAB",
            "GRACE",
            "GIST",
            "GROUP",
            "HBOC",
            "HELM",
            "HERSTORY",
            "HMHN",
            "HPT",
            "HTML",
            "HTTP",
            "ID",
            "INDEX",
            "INFO",
            "IREGXBSE",
            "JAVASCRIPT",
            "JIRA",
            "JPS",
            "JT",
            "JQUERY",
            "JSON",
            "JWT",
            "K2EDTA",
            "K3EDTA",
            "K8S",
            "LBJ",
            "LDAP",
            "LESS",
            "LOCAL",
            "LOG",
            "MALE",
            "MANS",
            "MAP",
            "METEOR",
            "METRIC",
            "MODE",
            "MONGO",
            "MYSQL",
            "NGS",
            "NODE",
            "NPM",
            "NRBA",
            "NREUM",
            "OAUTH",
            "OHSU",
            "OHSUMEL",
            "OSUPAN",
            "PAN",
            "PATCH",
            "PATIENT",
            "PDF",
            "PECGS",
            "POSTGRES",
            "POST",
            "PPAP",
            "PR",
            "PROD",
            "PUT",
            "QA",
            "REACT",
            "REDIS",
            "REMOTE",
            "REQ",
            "RESOURCE",
            "REST",
            "RESTORE",
            "SAML",
            "SASS",
            "SCSS",
            "SDK",
            "SLACK",
            "SMS",
            "SOAP",
            "SQL",
            "SSH",
            "SSL",
            "STAGE",
            "START",
            "TEAMS",
            "TERRAFORM",
            "TEST",
            "TLS",
            "TMDU",
            "TOKEN",
            "TRACE",
            "TYPESCRIPT",
            "UCLA",
            "UHN",
            "UI",
            "UID",
            "UNIFY",
            "URL",
            "USER",
            "UX",
            "VPN",
            "VUE",
            "WARN",
            "WEBPACK",
            "XHR",
            "XML",
            "YARN",
            "YASOJIMA",
            "YOSHIDA",
            "ZEPTO",
            "ZOOM",
        }

        if gene in web_terms:
            return False

        # Skip if it looks like a product/test code (starts with PR, contains only digits after letters)
        if gene.startswith("PR") and len(gene) > 2:
            return False

        # Skip very short codes that are likely not genes
        if len(gene) <= 2:
            return False

        # Skip if it's all numbers
        if gene.isdigit():
            return False

        # Must contain at least one letter that's not just at the beginning
        if len([c for c in gene if c.isalpha()]) < 2:
            return False

        return True

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
            chrome_options.add_argument("--remote-debugging-port=9222")
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

            # R script uses: //meta[contains(@content,"AIP")]
            # Look for meta tags containing "AIP" in their content
            meta_elements = driver.find_elements(By.TAG_NAME, "meta")
            for meta in meta_elements:
                content = meta.get_attribute("content") or ""
                if "AIP" in content:
                    import re

                    # Split content by commas and extract gene symbols (like R script)
                    potential_genes = content.split(",")
                    for gene in potential_genes:
                        gene = gene.strip()
                        if gene:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if (
                                cleaned_gene
                                and self.validate_gene_symbol(cleaned_gene)
                                and self._is_likely_gene(cleaned_gene)
                            ):
                                genes.append(cleaned_gene)
                    break

            # Look for genes in JavaScript variables with improved filtering
            if not genes:
                scripts = driver.find_elements(By.TAG_NAME, "script")
                for script in scripts:
                    script_content = script.get_attribute("innerHTML")
                    if script_content and any(
                        term in script_content for term in ["genes", "panel", "test"]
                    ):
                        import re

                        # Extract potential gene symbols with basic pattern
                        potential_genes = re.findall(
                            r"\b[A-Z][A-Z0-9]{1,7}\b", script_content
                        )

                        # Apply better filtering
                        for gene in potential_genes:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if (
                                cleaned_gene
                                and self.validate_gene_symbol(cleaned_gene)
                                and self._is_likely_gene(cleaned_gene)
                            ):
                                genes.append(cleaned_gene)

            # Look for genes in page content with better filtering
            if not genes:
                # Get the full page text and look for gene symbols
                page_text = driver.find_element(By.TAG_NAME, "body").text
                import re

                potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,7}\b", page_text)

                for gene in potential_genes:
                    cleaned_gene = self.clean_gene_symbol(gene)
                    if (
                        cleaned_gene
                        and self.validate_gene_symbol(cleaned_gene)
                        and self._is_likely_gene(cleaned_gene)
                    ):
                        genes.append(cleaned_gene)

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
            # Restore original no_proxy environment variable
            if "original_no_proxy" in locals():
                if original_no_proxy:
                    os.environ["no_proxy"] = original_no_proxy
                else:
                    os.environ.pop("no_proxy", None)
            if driver:
                driver.quit()
