"""
Eurogentest NGS panel SNP scraper.

Extracts identity SNPs from the Eurogentest NGS guidelines PDF.
Reference: https://www.college-genetics.be/assets/recommendations/fr/guidelines/EuroGentest%20NGS_2014.pdf
"""

import logging
import re
from pathlib import Path
from typing import Any

import pdfplumber

from custom_panel.sources_snp.downloader import PanelDownloader

from .base_snp_scraper import BaseSNPScraper

logger = logging.getLogger(__name__)


class EurogentestParser(BaseSNPScraper):
    """
    Parser for Eurogentest NGS identity SNP panel.

    This parser extracts rsIDs from the Eurogentest guidelines PDF that
    includes SNPs for sample tracking in NGS workflows.
    """

    def parse(self) -> dict[str, Any]:
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

            # Extract rsIDs and coordinates from PDF
            snp_data = self._extract_rsids_from_pdf(pdf_path)

            # Create SNP records with coordinates when available
            snps = []
            for data in snp_data:
                rsid = data["rsid"]

                # Extract all non-rsid fields for the record
                record_kwargs = {k: v for k, v in data.items() if k != "rsid"}

                snps.append(
                    self.create_snp_record(
                        rsid=rsid,
                        category="identity",
                        panel_specific_name="Eurogentest NGS 2014",
                        **record_kwargs,
                    ),
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from Eurogentest panel",
            )

            return {
                "panel_name": "eurogentest_ngs_panel",
                "source_url": self.url,
                "description": "Eurogentest NGS Guidelines (2014) - Identity SNP panel",
                "snps": snps,
                "metadata": {
                    "publication": "Eurogentest Guidelines 2014",
                    "organization": "European Society of Human Genetics",
                    "source_url": self.url,
                    "snp_count": len(snps),
                },
            }

        except Exception as e:
            logger.error(f"Error parsing Eurogentest panel: {e}")
            raise

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[dict[str, Any]]:
        """
        Extract rsIDs with chromosome information from the Eurogentest PDF.

        The PDF contains rsIDs in a table format with chromosome info.
        Pattern: chr.+rs[0-9]* (chromosome info followed by rsID)

        Args:
            pdf_path: Path to the PDF file

        Returns:
            List of dictionaries with rsID and coordinate information
        """
        snp_records = []

        try:
            with pdfplumber.open(pdf_path) as pdf:
                for page in pdf.pages:
                    text = page.extract_text()
                    if not text:
                        continue

                    # Split into lines
                    lines = text.split("\n")

                    for line in lines:
                        line = line.strip()

                        # Look for lines with chromosome info and rsID
                        # Pattern: chr followed by position info and rsID
                        if re.search(r"chr.+rs\d+", line, re.IGNORECASE):
                            # Extract both chromosome and rsID information
                            snp_record = self._parse_eurogentest_line(line)
                            if snp_record:
                                snp_records.append(snp_record)

                        # Also check for standalone rsIDs
                        elif re.match(r"^rs\d+", line, re.IGNORECASE):
                            rsid_match = re.match(r"^(rs\d+)", line, re.IGNORECASE)
                            if rsid_match:
                                rsid = self.clean_rsid(rsid_match.group(1))
                                if self.validate_rsid(rsid):
                                    snp_records.append({"rsid": rsid})

            # Remove duplicates based on rsID while preserving order
            seen = set()
            unique_records = []
            for record in snp_records:
                rsid = record["rsid"]
                if rsid not in seen:
                    seen.add(rsid)
                    unique_records.append(record)

            with_coords = sum(1 for r in unique_records if "chromosome" in r)
            logger.info(
                f"Extracted {len(unique_records)} SNPs ({with_coords} with chromosome info)",
            )
            return unique_records

        except Exception as e:
            logger.error(f"Error extracting rsIDs from PDF: {e}")
            raise

    def _parse_eurogentest_line(self, line: str) -> dict[str, Any] | None:
        """
        Parse a line from Eurogentest PDF with chromosome and rsID information.

        Expected pattern: chr followed by position info and rsID
        Examples:
        - "chr1 123456789 rs1234567"
        - "1:123456789 rs1234567"

        Args:
            line: Text line from PDF

        Returns:
            Dictionary with parsed SNP information, or None if parsing fails
        """
        try:
            # Extract rsID first
            rsid_match = re.search(r"\b(rs\d+)\b", line)
            if not rsid_match:
                return None

            rsid = self.clean_rsid(rsid_match.group(1))
            if not self.validate_rsid(rsid):
                return None

            record = {"rsid": rsid}

            # Extract chromosome information
            # Look for chr pattern followed by chromosome number/letter
            chr_match = re.search(r"chr\s*([0-9XYM]+)", line, re.IGNORECASE)
            if chr_match:
                chromosome = chr_match.group(1)
                record["chromosome"] = self.clean_chromosome(chromosome)
            else:
                # Look for standalone chromosome notation
                chr_match = re.search(r"\b([0-9]+|X|Y|M|MT)\s*:", line, re.IGNORECASE)
                if chr_match:
                    chromosome = chr_match.group(1)
                    record["chromosome"] = self.clean_chromosome(chromosome)

            return record

        except Exception as e:
            logger.debug(f"Failed to parse line: {e}")
            return None
