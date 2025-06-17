"""
Fulgent Genetics parser.

This module extracts gene information from Fulgent Genetics panel pages.
Now uses BeautifulSoup to parse the content-details section containing gene lists.
"""

import logging
import os
import re

import requests
from bs4 import BeautifulSoup

from .base_parser import BaseParser

logger = logging.getLogger(__name__)


class FulgentParser(BaseParser):
    """Parser for Fulgent Genetics panels."""

    def parse(self) -> list[str]:
        """
        Parse Fulgent Genetics page and extract gene symbols.

        Returns:
            List of gene symbols found on the page

        Raises:
            Exception: If parsing fails
        """
        try:
            # Setup proxy configuration for Charite network
            proxies = {}
            if os.environ.get("http_proxy") or os.environ.get("https_proxy"):
                proxies = {
                    "http": os.environ.get("http_proxy", "http://proxy.charite.de:8080"),
                    "https": os.environ.get("https_proxy", "http://proxy.charite.de:8080")
                }
                logger.info(f"Using proxy configuration: {proxies}")
            
            # Setup headers to mimic a regular browser
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
                'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.9',
                'Accept-Encoding': 'gzip, deflate, br',
                'DNT': '1',
                'Connection': 'keep-alive',
                'Upgrade-Insecure-Requests': '1',
            }
            
            # Make request with proxy support and browser headers
            try:
                response = requests.get(
                    self.url, 
                    timeout=30,
                    headers=headers,
                    proxies=proxies if proxies else None
                )
                response.raise_for_status()
            except requests.RequestException as e:
                logger.error(
                    f"Network request to Fulgent Genetics URL failed: {self.url} - {e}"
                )
                raise

            # Parse HTML
            soup = BeautifulSoup(response.content, "html.parser")
            genes = []

            # Look for the content-details section containing genes
            content_details = soup.find("div", class_="content-details")
            if content_details:
                logger.info("Found content-details section")
                
                # Look for text that contains "Genes:" followed by comma-separated gene list
                text = content_details.get_text()
                
                # Find the genes section after "Genes:" label
                # Look for both "Genes:" and "Genes analyzed:" patterns
                genes_match = re.search(r"Genes\s*(?:analyzed)?:\s*([A-Z0-9,\s]+)", text, re.IGNORECASE | re.MULTILINE)
                if genes_match:
                    genes_text = genes_match.group(1)
                    logger.info(f"Found genes text: {genes_text[:100]}...")
                    
                    # Parse comma-separated genes
                    gene_parts = genes_text.split(",")
                    for gene_part in gene_parts:
                        gene_part = gene_part.strip()
                        cleaned_gene = self.clean_gene_symbol(gene_part)
                        if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                            genes.append(cleaned_gene)
                else:
                    # Fallback: Look for long comma-separated lists of gene-like strings
                    logger.info("Genes: label not found, trying fallback approach")
                    
                    # Look for patterns like "AIP, ALK, APC, ATM, AXIN2, BAP1, BARD1, BMPR1A, BRCA1, BRCA2"
                    # Find sequences of comma-separated uppercase words
                    comma_separated_match = re.search(r"\b([A-Z][A-Z0-9]{1,7}(?:,\s*[A-Z][A-Z0-9]{1,7}){10,})", text)
                    if comma_separated_match:
                        genes_text = comma_separated_match.group(1)
                        logger.info(f"Found comma-separated genes: {genes_text[:100]}...")
                        
                        gene_parts = genes_text.split(",")
                        for gene_part in gene_parts:
                            gene_part = gene_part.strip()
                            cleaned_gene = self.clean_gene_symbol(gene_part)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                    else:
                        # Extract all potential gene symbols but be more restrictive
                        potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{2,7}\b", text)
                        # Filter out common false positives
                        filtered_genes = []
                        false_positives = {"RISE", "STAT", "VUS", "CGI", "FULL", "TEST", "DNA", "RNA", "PCR", "NGS", "FDA", "USA", "LLC", "INC", "CORP", "LABS", "LAB", "TECH", "INFO", "CARE", "HEALTH", "MED", "MEDICAL", "CLINIC", "CENTER", "GENE", "GENES", "PANEL", "PANELS", "COMPREHENSIVE", "CANCER", "TUMOR", "TUMORS"}
                        
                        for gene in potential_genes:
                            if gene not in false_positives:
                                cleaned_gene = self.clean_gene_symbol(gene)
                                if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                    filtered_genes.append(cleaned_gene)
                        
                        # Only use this approach if we get a reasonable number of genes
                        if len(filtered_genes) > 10:
                            genes.extend(filtered_genes)
                            logger.info(f"Found {len(filtered_genes)} genes using fallback approach")

            # If content-details section not found, try other approaches
            if not genes:
                logger.info("Content-details section not found, trying meta tags")
                
                # Look in meta tags (original approach)
                meta_tags = soup.find_all("meta")
                for meta in meta_tags:
                    content = meta.get("content", "")
                    if "AIP" in content or "gene" in content.lower():
                        # Parse genes from content
                        potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,7}\b", content)
                        for gene in potential_genes:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)

            # Final fallback: look in any element containing gene-related keywords
            if not genes:
                logger.info("Meta tags approach failed, trying full page scan")
                for element in soup.find_all(["div", "p", "span", "section"]):
                    text = element.get_text()
                    if any(
                        keyword in text.lower()
                        for keyword in ["genes:", "gene list", "panel genes"]
                    ):
                        potential_genes = re.findall(r"\b[A-Z][A-Z0-9]{1,7}\b", text)
                        for gene in potential_genes:
                            cleaned_gene = self.clean_gene_symbol(gene)
                            if cleaned_gene and self.validate_gene_symbol(cleaned_gene):
                                genes.append(cleaned_gene)
                        if genes:  # Stop after first successful match
                            break

            # Remove duplicates while preserving order
            unique_genes = []
            seen = set()
            for gene in genes:
                if gene not in seen:
                    unique_genes.append(gene)
                    seen.add(gene)

            logger.info(
                f"Extracted {len(unique_genes)} genes from Fulgent Genetics panel"
            )
            return unique_genes

        except Exception as e:
            logger.error(f"Error parsing Fulgent Genetics panel: {e}")
            raise
