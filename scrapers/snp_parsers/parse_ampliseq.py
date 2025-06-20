"""
Ampliseq Illumina Sample ID panel SNP scraper.

Extracts identity SNPs from the Illumina Ampliseq Sample ID panel PDF.
Reference: https://support.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/ampliseq-sample-id-panel-app-note-m-gl-01513/ampliseq-sample-id-panel-app-note-m-gl-01513.pdf
"""

import logging
import re
from pathlib import Path
from typing import Any

import pdfplumber

from custom_panel.sources_snp.downloader import PanelDownloader

from .base_snp_scraper import BaseSNPScraper

logger = logging.getLogger(__name__)


class AmpliseqParser(BaseSNPScraper):
    """
    Parser for Illumina Ampliseq Sample ID panel.

    This parser extracts rsIDs from the Illumina PDF that describes
    their amplicon-based sample identification panel.
    """

    def parse(self) -> dict[str, Any]:
        """
        Parse the Ampliseq PDF and extract SNP data.

        Returns:
            Dictionary containing panel data and SNPs
        """
        logger.info(f"Parsing Ampliseq panel from {self.url}")

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
                        panel_specific_name="Illumina AmpliSeq Sample ID Panel",
                    )
                )

            logger.info(f"Successfully extracted {len(snps)} SNPs from Ampliseq panel")

            return {
                "panel_name": "ampliseq_illumina_panel",
                "source_url": self.url,
                "description": "Illumina AmpliSeq Sample ID Panel - Amplicon-based sample tracking",
                "snps": snps,
                "metadata": {
                    "vendor": "Illumina",
                    "technology": "AmpliSeq amplicon sequencing",
                    "source_url": self.url,
                    "snp_count": len(snps),
                },
            }

        except Exception as e:
            logger.error(f"Error parsing Ampliseq panel: {e}")
            raise

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[str]:
        """
        Extract rsIDs from the Ampliseq PDF.

        The PDF contains rsIDs in a table format.
        Pattern: ^[0-9].+rs[0-9]* (number at start, rsID somewhere in line)

        Args:
            pdf_path: Path to the PDF file

        Returns:
            List of unique rsIDs
        """
        rsids = []

        try:
            with pdfplumber.open(pdf_path) as pdf:
                for page in pdf.pages:
                    text = page.extract_text()
                    if not text:
                        continue

                    # Try to extract tables first
                    tables = page.extract_tables()
                    for table in tables:
                        for row in table:
                            if row:
                                # Check each cell in the row
                                for cell in row:
                                    if cell:
                                        rsids.extend(self.extract_rsids_from_text(cell))

                    # Also check plain text
                    lines = text.split("\n")

                    for line in lines:
                        line = line.strip()

                        # Look for lines starting with number and containing rsID
                        # Pattern: number at start, rsID somewhere in line
                        if re.match(r"^\d+.+rs\d+", line):
                            # Extract rsID from line
                            rsid_match = re.search(r"(rs\d+)", line)
                            if rsid_match:
                                rsid = rsid_match.group(1)
                                if self.validate_rsid(rsid):
                                    rsids.append(rsid)

                        # Also check for standalone rsIDs
                        elif re.match(r"^rs\d+", line):
                            rsid_match = re.match(r"^(rs\d+)", line)
                            if rsid_match:
                                rsid = rsid_match.group(1)
                                if self.validate_rsid(rsid):
                                    rsids.append(rsid)

                        # Extract any rsIDs from the line
                        else:
                            line_rsids = self.extract_rsids_from_text(line)
                            rsids.extend(line_rsids)

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
