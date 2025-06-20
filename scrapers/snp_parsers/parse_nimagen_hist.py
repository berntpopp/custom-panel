"""
NimaGen Human Identification and Sample Tracking (HIST) kit SNP scraper.

Extracts identity SNPs from the NimaGen HIST kit PDF.
Reference: https://www.nimagen.com/gfx/brands/NGS%20Flyers%20and%20Folders/HEST%20Brochure.pdf
"""

import logging
import re
from pathlib import Path
from typing import Any

import pdfplumber

from custom_panel.sources_snp.downloader import PanelDownloader

from .base_snp_scraper import BaseSNPScraper

logger = logging.getLogger(__name__)


class NimaGenHISTParser(BaseSNPScraper):
    """
    Parser for NimaGen Human Identification and Sample Tracking kit.

    This parser extracts rsIDs from the PDF that contains a table with
    columns: # dbSNP Chr# Location HG38 Gene MAF 1000 genomes
    with 34 SNPs for sample tracking.
    """

    def parse(self) -> dict[str, Any]:
        """
        Parse the NimaGen HIST PDF and extract SNP data.

        Returns:
            Dictionary containing panel data and SNPs
        """
        logger.info(f"Parsing NimaGen HIST panel from {self.url}")

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
                        panel_specific_name="NimaGen Human Identification and Sample Tracking Kit",
                        **record_kwargs,
                    )
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from NimaGen HIST panel"
            )

            return {
                "panel_name": "nimagen_hist_panel",
                "source_url": self.url,
                "description": "NimaGen Human Identification and Sample Tracking Kit - 34 SNP panel",
                "snps": snps,
                "metadata": {
                    "vendor": "NimaGen",
                    "kit_name": "Human Identification and Sample Tracking Kit",
                    "technology": "NGS-based sample tracking",
                    "source_url": self.url,
                    "snp_count": len(snps),
                },
            }

        except Exception as e:
            logger.error(f"Error parsing NimaGen HIST panel: {e}")
            raise

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[dict[str, Any]]:
        """
        Extract rsIDs with genomic coordinates from the NimaGen HIST PDF.

        The PDF contains a table with columns:
        # dbSNP Chr# Location HG38 Gene MAF 1000 genomes

        Args:
            pdf_path: Path to the PDF file

        Returns:
            List of dictionaries with rsID and coordinate information
        """
        snp_records = []
        target_table_found = False

        try:
            with pdfplumber.open(pdf_path) as pdf:
                for page_num, page in enumerate(pdf.pages):
                    # Check if this page contains SNP data
                    text = page.extract_text()
                    if not text:
                        continue

                    # Look for table headers indicating SNP data
                    if (
                        (
                            "dbSNP" in text
                            and "Chr#" in text
                            and "Location" in text
                            and "HG38" in text
                        )
                        or ("# dbSNP Chr# Location HG38" in text)
                        or (
                            "Human Identification" in text and "Sample Tracking" in text
                        )
                    ):
                        logger.info(f"Found potential SNP table on page {page_num + 1}")
                        target_table_found = True

                        # Try to extract tables from this page
                        tables = page.extract_tables()
                        for table in tables:
                            if not table:
                                continue

                            # Look for table with SNP data
                            # Expected columns: # dbSNP Chr# Location HG38 Gene MAF 1000 genomes
                            for row_idx, row in enumerate(table):
                                if not row:
                                    continue

                                # Skip header rows - look for rows starting with numbers
                                if row_idx == 0:
                                    # Check if this looks like our target header
                                    header_text = " ".join(
                                        str(cell) for cell in row if cell
                                    )
                                    if (
                                        "dbSNP" in header_text
                                        and "Chr#" in header_text
                                        and "Location" in header_text
                                    ):
                                        logger.info(
                                            f"Found HIST coordinate table header: {header_text}"
                                        )
                                    continue

                                # Check if first column looks like a row number
                                if not row[0]:
                                    continue

                                try:
                                    int(str(row[0]).strip())
                                except (ValueError, AttributeError):
                                    continue

                                # Parse the full row for genomic information
                                snp_record = self._parse_hist_table_row(row)
                                if snp_record:
                                    snp_records.append(snp_record)
                                else:
                                    # Debug failed parsing
                                    row_text = [
                                        str(cell) for cell in row[:6]
                                    ]  # First 6 columns
                                    logger.debug(f"Failed to parse row: {row_text}")

                        # Also try text extraction for this page
                        lines = text.split("\n")
                        for line in lines:
                            # Look for lines that contain rsIDs
                            rsid_matches = re.findall(r"\b(rs\d+)\b", line)
                            for rsid in rsid_matches:
                                if self.validate_rsid(rsid):
                                    snp_records.append({"rsid": rsid})

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
                f"Extracted {len(unique_records)} SNPs ({with_coords} with coordinates) from NimaGen HIST PDF"
            )
            return unique_records

        except Exception as e:
            logger.error(f"Error extracting rsIDs from PDF: {e}")
            raise

    def _parse_hist_table_row(self, row: list) -> dict[str, Any] | None:
        """
        Parse a single row from the NimaGen HIST table.

        Expected format: # dbSNP Chr# Location HG38 Gene MAF 1000 genomes
        Example: 1 rs1410592 1 179551371 NPHS2 , AXDND1 G 0.4123

        Args:
            row: Table row as list of cells

        Returns:
            Dictionary with parsed SNP information, or None if parsing fails
        """
        try:
            if len(row) < 6:  # Need at least 6 columns for basic info
                return None

            # Extract data from expected positions
            rsid = str(row[1]).strip() if row[1] else ""
            chromosome = str(row[2]).strip() if row[2] else ""
            position_hg38 = str(row[3]).strip() if row[3] else ""
            gene = str(row[4]).strip() if row[4] else ""

            # Validate rsID
            if not self.validate_rsid(rsid):
                return None

            record = {"rsid": rsid}

            # Add chromosome if available
            if chromosome:
                record["chromosome"] = self.clean_chromosome(chromosome)

            # Parse HG38 position (plain integer format)
            if position_hg38 and position_hg38.isdigit():
                try:
                    position = int(position_hg38)
                    record["position"] = position
                    record["assembly"] = "GRCh38"
                except ValueError:
                    logger.debug(f"Could not parse position: {position_hg38}")

            # Add gene information (handle comma-separated genes)
            if gene:
                # Clean up gene names - remove extra spaces and handle comma separation
                genes = [g.strip() for g in gene.split(",") if g.strip()]
                if genes:
                    record["gene"] = genes[0] if len(genes) == 1 else genes

            return record

        except Exception as e:
            logger.debug(f"Failed to parse table row: {e}")
            return None
