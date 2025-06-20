"""
NimaGen Human DNA Sample Identification Kit (HEST) SNP scraper.

Extracts identity SNPs from the NimaGen HEST kit PDF.
Reference: https://www.nimagen.com/gfx/HEST/NimaGen_folder%20EasySeq%20Human%20DNA%20Sample%20Identification%20Kits%20v1-1.pdf
"""

import logging
import re
from pathlib import Path
from typing import Any

import pdfplumber

from custom_panel.sources_snp.downloader import PanelDownloader

from .base_snp_scraper import BaseSNPScraper

logger = logging.getLogger(__name__)


class NimaGenHESTParser(BaseSNPScraper):
    """
    Parser for NimaGen EasySeq Human DNA Sample Identification Kit.

    This parser extracts rsIDs from the PDF that contains
    "Table 1 | EasySeq™ Human DNA Sample Identification Kit Targets Overview"
    with 37 SNPs for sample identification.
    """

    def parse(self) -> dict[str, Any]:
        """
        Parse the NimaGen HEST PDF and extract SNP data.

        Returns:
            Dictionary containing panel data and SNPs
        """
        logger.info(f"Parsing NimaGen HEST panel from {self.url}")

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
                        panel_specific_name="NimaGen EasySeq Human DNA Sample Identification Kit",
                    )
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from NimaGen HEST panel"
            )

            return {
                "panel_name": "nimagen_hest_panel",
                "source_url": self.url,
                "description": "NimaGen EasySeq Human DNA Sample Identification Kit - 37 SNP panel",
                "snps": snps,
                "metadata": {
                    "vendor": "NimaGen",
                    "kit_name": "EasySeq Human DNA Sample Identification Kit",
                    "technology": "NGS-based sample identification",
                    "source_url": self.url,
                    "snp_count": len(snps),
                },
            }

        except Exception as e:
            logger.error(f"Error parsing NimaGen HEST panel: {e}")
            raise

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[str]:
        """
        Extract rsIDs from the NimaGen HEST PDF.

        The PDF contains "Table 1 | EasySeq™ Human DNA Sample Identification Kit Targets Overview"
        with columns: # SNP Chr Gene Location HG38 Location HG19 MAF ALFA Total

        Args:
            pdf_path: Path to the PDF file

        Returns:
            List of unique rsIDs
        """
        rsids = []
        target_table_found = False

        try:
            with pdfplumber.open(pdf_path) as pdf:
                for page_num, page in enumerate(pdf.pages):
                    # Check if this page contains the target table
                    text = page.extract_text()
                    if not text:
                        continue

                    # Look for the specific table title
                    if (
                        "EasySeq" in text
                        and "Human DNA Sample Identification" in text
                        and "Targets Overview" in text
                    ):
                        logger.info(f"Found target table on page {page_num + 1}")
                        target_table_found = True

                        # Try to extract tables from this page
                        tables = page.extract_tables()
                        for table in tables:
                            if not table:
                                continue

                            # Look for table with SNP data
                            for row_idx, row in enumerate(table):
                                if not row:
                                    continue

                                # Skip header rows - look for rows starting with numbers
                                if row_idx == 0 or not row[0]:
                                    continue

                                # Check if first column looks like a row number
                                try:
                                    int(str(row[0]).strip())
                                except (ValueError, AttributeError):
                                    continue

                                # Look for rsID in the row (typically second column)
                                for cell in row[
                                    1:4
                                ]:  # Check first few columns after row number
                                    if cell and isinstance(cell, str):
                                        rsid_match = re.search(r"\b(rs\d+)\b", cell)
                                        if rsid_match:
                                            rsid = rsid_match.group(1)
                                            if self.validate_rsid(rsid):
                                                rsids.append(rsid)
                                                break

                        # Also try text extraction for this page
                        lines = text.split("\n")
                        for line in lines:
                            # Look for lines that contain rsIDs
                            rsid_matches = re.findall(r"\b(rs\d+)\b", line)
                            for rsid in rsid_matches:
                                if self.validate_rsid(rsid):
                                    rsids.append(rsid)

            # If we didn't find the specific table, try extracting from entire document
            if not target_table_found:
                logger.warning(
                    "Target table not found, extracting rsIDs from entire document"
                )
                with pdfplumber.open(pdf_path) as pdf:
                    for page in pdf.pages:
                        text = page.extract_text()
                        if text:
                            rsid_matches = re.findall(r"\b(rs\d+)\b", text)
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

            logger.info(
                f"Extracted {len(unique_rsids)} unique rsIDs from NimaGen HEST PDF"
            )
            return unique_rsids

        except Exception as e:
            logger.error(f"Error extracting rsIDs from PDF: {e}")
            raise
