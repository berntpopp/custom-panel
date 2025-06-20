"""
Eurogentest NGS panel SNP scraper.

Extracts identity SNPs from the Eurogentest NGS guidelines PDF.
Reference: https://www.college-genetics.be/assets/recommendations/fr/guidelines/EuroGentest%20NGS_2014.pdf
"""

import logging
from typing import Any, Dict, List
import pdfplumber
import re
from pathlib import Path

from .base_snp_scraper import BaseSNPScraper
from custom_panel.sources_snp.downloader import PanelDownloader

logger = logging.getLogger(__name__)


class EurogentestParser(BaseSNPScraper):
    """
    Parser for Eurogentest NGS identity SNP panel.
    
    This parser extracts rsIDs from the Eurogentest guidelines PDF that
    includes SNPs for sample tracking in NGS workflows.
    """

    def parse(self) -> Dict[str, Any]:
        """
        Parse the Eurogentest PDF and extract SNP data.

        Returns:
            Dictionary containing panel data and SNPs
        """
        logger.info(f"Parsing Eurogentest panel from {self.url}")
        
        try:
            # Download PDF
            with PanelDownloader() as downloader:
                pdf_path = downloader.download(self.url, "pdf", "requests")
            
            # Extract rsIDs from PDF
            rsids = self._extract_rsids_from_pdf(pdf_path)
            
            # Create SNP records
            snps = []
            for rsid in rsids:
                snps.append(self.create_snp_record(
                    rsid=rsid,
                    category="identity",
                    panel_specific_name="Eurogentest NGS 2014"
                ))
            
            logger.info(f"Successfully extracted {len(snps)} SNPs from Eurogentest panel")
            
            return {
                "panel_name": "eurogentest_ngs_panel",
                "source_url": self.url,
                "description": "Eurogentest NGS Guidelines (2014) - Identity SNP panel",
                "snps": snps,
                "metadata": {
                    "publication": "Eurogentest Guidelines 2014",
                    "organization": "European Society of Human Genetics",
                    "source_url": self.url,
                    "snp_count": len(snps)
                }
            }
            
        except Exception as e:
            logger.error(f"Error parsing Eurogentest panel: {e}")
            raise

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> List[str]:
        """
        Extract rsIDs from the Eurogentest PDF.
        
        The PDF contains rsIDs in a table format with chromosome info.
        Pattern: chr.+rs[0-9]* (chromosome info followed by rsID)
        
        Args:
            pdf_path: Path to the PDF file
            
        Returns:
            List of unique rsIDs
        """
        rsids = []
        
        try:
            with pdfplumber.open(pdf_path) as pdf:
                for page_num, page in enumerate(pdf.pages):
                    text = page.extract_text()
                    if not text:
                        continue
                    
                    # Split into lines
                    lines = text.split('\n')
                    
                    for line in lines:
                        line = line.strip()
                        
                        # Look for lines with chromosome info and rsID
                        # Pattern: chr followed by position info and rsID
                        if re.search(r'chr.+rs\d+', line, re.IGNORECASE):
                            # Extract rsID from line
                            # Split by whitespace and look for rsID pattern
                            parts = line.split()
                            for part in parts:
                                if re.match(r'^rs\d+$', part, re.IGNORECASE):
                                    rsid = self.clean_rsid(part)
                                    if self.validate_rsid(rsid):
                                        rsids.append(rsid)
                        
                        # Also check for standalone rsIDs
                        elif re.match(r'^rs\d+', line, re.IGNORECASE):
                            rsid_match = re.match(r'^(rs\d+)', line, re.IGNORECASE)
                            if rsid_match:
                                rsid = self.clean_rsid(rsid_match.group(1))
                                if self.validate_rsid(rsid):
                                    rsids.append(rsid)
            
            # Remove duplicates while preserving order
            seen = set()
            unique_rsids = []
            for rsid in rsids:
                if rsid not in seen:
                    seen.add(rsid)
                    unique_rsids.append(rsid)
            
            return unique_rsids
            
        except Exception as e:
            logger.error(f"Error extracting rsIDs from PDF: {e}")
            raise