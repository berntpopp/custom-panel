"""
Mayo Clinic parser.

This module extracts gene information from Mayo Clinic panel pages.
Looks for genes in italic tags.
"""

import logging

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class MayoParser(BaseParser):
    """Parser for Mayo Clinic panels."""

    def parse(self) -> list[str]:
        """
        Parse Mayo Clinic panel page and extract gene symbols.

        Returns:
            List of gene symbols found on the page

        Raises:
            Exception: If parsing fails
        """
        try:
            # Try requests first, fallback to Selenium for strong anti-bot protection
            html_content = None
            try:
                import time
                
                headers = {
                    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
                    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
                    'Accept-Language': 'en-US,en;q=0.9',
                    'Accept-Encoding': 'gzip, deflate, br',
                    'DNT': '1',
                    'Connection': 'keep-alive',
                    'Upgrade-Insecure-Requests': '1',
                    'Sec-Fetch-Dest': 'document',
                    'Sec-Fetch-Mode': 'navigate',
                    'Sec-Fetch-Site': 'cross-site',
                    'Sec-Fetch-User': '?1',
                    'Cache-Control': 'max-age=0',
                    'sec-ch-ua': '"Not_A Brand";v="8", "Chromium";v="120", "Google Chrome";v="120"',
                    'sec-ch-ua-mobile': '?0',
                    'sec-ch-ua-platform': '"macOS"'
                }
                
                session = requests.Session()
                session.headers.update(headers)
                response = session.get(self.url, timeout=30)
                response.raise_for_status()
                html_content = response.content
                
            except requests.RequestException as e:
                logger.warning(f"Requests method failed for Mayo Clinic: {e}. Trying Selenium...")
                
                # Fallback to Selenium
                try:
                    from selenium import webdriver
                    from selenium.webdriver.chrome.options import Options
                    from selenium.webdriver.common.by import By
                    from selenium.webdriver.support.ui import WebDriverWait
                    from selenium.webdriver.support import expected_conditions as EC
                    import time
                    
                    # Setup Chrome options
                    chrome_options = Options()
                    chrome_options.add_argument('--headless')  # Run in background
                    chrome_options.add_argument('--no-sandbox')
                    chrome_options.add_argument('--disable-dev-shm-usage')
                    chrome_options.add_argument('--disable-gpu')
                    chrome_options.add_argument('--window-size=1920,1080')
                    chrome_options.add_argument('--user-agent=Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36')
                    chrome_options.add_argument('--no-proxy-server')  # Bypass proxy for direct connection
                    chrome_options.add_argument('--disable-web-security')
                    chrome_options.add_argument('--allow-running-insecure-content')
                    chrome_options.add_argument('--disable-features=VizDisplayCompositor')
                    
                    # Create driver
                    driver = webdriver.Chrome(options=chrome_options)
                    
                    try:
                        # Navigate to page
                        driver.get(self.url)
                        
                        # Wait for page to load
                        WebDriverWait(driver, 10).until(
                            EC.presence_of_element_located((By.TAG_NAME, "body"))
                        )
                        
                        # Additional wait for content to load
                        time.sleep(3)
                        
                        html_content = driver.page_source.encode('utf-8')
                        logger.info("Successfully retrieved Mayo Clinic page with Selenium")
                        
                    finally:
                        driver.quit()
                        
                except Exception as selenium_error:
                    logger.error(f"Selenium method also failed for Mayo Clinic: {selenium_error}")
                    raise Exception(f"Both requests and Selenium failed. Last error: {selenium_error}")
            
            if not html_content:
                raise Exception("Failed to retrieve page content")

            # Parse HTML
            soup = BeautifulSoup(html_content, "html.parser")

            genes = []

            # Look for genes in the Genetics Test Information section
            genetics_info = soup.find("div", class_="GeneticsTestInformation-value")
            if genetics_info:
                # Find the paragraph containing the gene list
                paragraphs = genetics_info.find_all("p")
                for paragraph in paragraphs:
                    text = paragraph.get_text(strip=True)
                    # Check if this paragraph contains gene symbols (comma-separated list)
                    if "," in text and any(
                        word.isupper() and len(word) <= 7 for word in text.split()
                    ):
                        logger.info(
                            "Found gene list in Genetics Test Information section"
                        )
                        # Split by comma and extract genes
                        gene_parts = text.split(",")
                        for gene_part in gene_parts:
                            # Clean up annotations like "(including promoters 1A and 1B)"
                            gene_part = gene_part.strip()
                            # Remove parenthetical information
                            import re

                            gene_part = re.sub(r"\s*\([^)]*\)", "", gene_part)

                            cleaned_gene = self.clean_gene_symbol(gene_part)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                        break  # Found the gene list, no need to continue

            # Fallback: Find italic tags containing gene lists (original logic)
            if not genes:
                italic_tags = soup.find_all("i")
                for italic in italic_tags:
                    text = italic.get_text(strip=True)
                    if text and (
                        "," in text or len(text.split()) > 1
                    ):  # Likely a gene list
                        # Split by comma and extract genes
                        for gene_part in text.split(","):
                            cleaned_gene = self.clean_gene_symbol(gene_part)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # Also check em tags
            if not genes:
                em_tags = soup.find_all("em")
                for em in em_tags:
                    text = em.get_text(strip=True)
                    if text and ("," in text or len(text.split()) > 1):
                        for gene_part in text.split(","):
                            cleaned_gene = self.clean_gene_symbol(gene_part)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # General search for gene patterns in divs/spans that might contain gene lists
            if not genes:
                for element in soup.find_all(["div", "span", "p", "li"]):
                    text = element.get_text()
                    # Look for comma-separated lists that might be genes
                    if "," in text and len(text) < 2000:  # Reasonable length limit
                        import re

                        # Extract potential gene symbols (2-7 uppercase letters/numbers)
                        potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,6}\b", text)
                        # Only proceed if we find multiple potential genes
                        if len(potential_genes) >= 5:
                            for gene in potential_genes:
                                cleaned_gene = self.clean_gene_symbol(gene)
                                if cleaned_gene and self.validate_gene_symbol(
                                    cleaned_gene
                                ):
                                    genes.append(cleaned_gene)

            # Remove duplicates while preserving order
            unique_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)

            logger.info(f"Extracted {len(unique_genes)} genes from Mayo Clinic panel")
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Mayo Clinic panel: {e}")
            raise
