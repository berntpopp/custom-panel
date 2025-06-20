"""
Pengelly panel SNP scraper.

Extracts identity SNPs from the Pengelly et al. paper PDF.
Reference: https://genomemedicine.biomedcentral.com/counter/pdf/10.1186/gm492.pdf
"""

import logging
import re
from pathlib import Path
from typing import Any

import pdfplumber

from custom_panel.sources_snp.downloader import PanelDownloader

from .base_snp_scraper import BaseSNPScraper

logger = logging.getLogger(__name__)


class PengellyParser(BaseSNPScraper):
    """
    Parser for Pengelly et al. identity SNP panel.

    This parser extracts rsIDs from the PDF publication that describes
    a panel of 56 SNPs for sample tracking in NGS studies.
    """

    def parse(self) -> dict[str, Any]:
        """
        Parse the Pengelly PDF and extract SNP data.

        Returns:
            Dictionary containing panel data and SNPs
        """
        logger.info(f"Parsing Pengelly panel from {self.url}")

        try:
            # Download PDF
            with PanelDownloader() as downloader:
                pdf_path = downloader.download(self.url, "pdf", "requests")

            # Extract rsIDs from PDF
            rsids = self._extract_rsids_from_pdf(pdf_path)

            # Create SNP records
            snps = []
            for rsid in rsids:
                snps.append(
                    self.create_snp_record(
                        rsid=rsid,
                        category="identity",
                        panel_specific_name="Pengelly et al. 2013",
                    )
                )

            logger.info(f"Successfully extracted {len(snps)} SNPs from Pengelly panel")

            return {
                "panel_name": "pengelly_panel",
                "source_url": self.url,
                "description": "Pengelly et al. (2013) - 56 SNP panel for sample tracking",
                "snps": snps,
                "metadata": {
                    "publication": "Genome Medicine 2013",
                    "doi": "10.1186/gm492",
                    "source_url": self.url,
                    "snp_count": len(snps),
                },
            }

        except Exception as e:
            logger.error(f"Error parsing Pengelly panel: {e}")
            raise

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[str]:
        """
        Extract rsIDs from the Pengelly PDF.

        The PDF contains rsIDs in table format, not at the start of lines.
        We need to extract all rsID patterns from the entire text.

        Args:
            pdf_path: Path to the PDF file

        Returns:
            List of unique rsIDs
        """
        rsids = []

        try:
            with pdfplumber.open(pdf_path) as pdf:
                # First try to extract tables
                for page in pdf.pages:
                    # Try table extraction first
                    tables = page.extract_tables()
                    for table in tables:
                        for row in table:
                            if row:
                                for cell in row:
                                    if cell and isinstance(cell, str):
                                        # Look for rsID in cell
                                        rsid_match = re.search(r"\b(rs\d+)\b", cell)
                                        if rsid_match:
                                            rsid = rsid_match.group(1)
                                            if self.validate_rsid(rsid):
                                                rsids.append(rsid)

                    # Also extract from plain text
                    text = page.extract_text()
                    if text:
                        # Extract all rsIDs from text
                        rsid_matches = re.findall(r"\b(rs\d+)\b", text, re.IGNORECASE)
                        for rsid in rsid_matches:
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
