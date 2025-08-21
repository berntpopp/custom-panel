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
                        panel_specific_name="NimaGen EasySeq Human DNA Sample Identification Kit",
                        **record_kwargs,
                    ),
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from NimaGen HEST panel",
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

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[dict[str, Any]]:
        """
        Extract rsIDs with genomic coordinates from the NimaGen HEST PDF.

        The PDF contains "Table 1 | EasySeq™ Human DNA Sample Identification Kit Targets Overview"
        with columns: # SNP Chr Gene Location HG38 Location HG19 MAF ALFA Total

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
                            # Expected columns: # SNP Chr Gene Location HG38 Location HG19 MAF ALFA Total
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

                                # Parse the full row for genomic information
                                snp_record = self._parse_hest_table_row(row)
                                if snp_record:
                                    snp_records.append(snp_record)

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
                    "Target table not found, extracting rsIDs from entire document",
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
                f"Extracted {len(unique_records)} SNPs ({with_coords} with coordinates) from NimaGen HEST PDF",
            )
            return unique_records

        except Exception as e:
            logger.error(f"Error extracting rsIDs from PDF: {e}")
            raise

    def _parse_hest_table_row(self, row: list) -> dict[str, Any] | None:
        """
        Parse a single row from the NimaGen HEST table.

        Expected format: # SNP Chr Gene Location HG38 Location HG19 MAF ALFA Total
        Example: 1 rs10203363 1 PTPN22 1:114377568 1:114377568 0.13 0.13 ...

        Args:
            row: Table row as list of cells

        Returns:
            Dictionary with parsed SNP information, or None if parsing fails
        """
        try:
            if len(row) < 7:  # Need at least 7 columns for basic info
                return None

            # Extract data from expected positions
            rsid = str(row[1]).strip() if row[1] else ""
            chromosome = str(row[2]).strip() if row[2] else ""
            gene = str(row[3]).strip() if row[3] else ""
            location_hg38 = str(row[4]).strip() if row[4] else ""
            location_hg19 = str(row[5]).strip() if row[5] else ""

            # Validate rsID
            if not self.validate_rsid(rsid):
                return None

            record = {"rsid": rsid}

            # Add chromosome if available
            if chromosome:
                record["chromosome"] = chromosome

            # Parse HG38 coordinates
            if location_hg38:
                hg38_coords = self.parse_genomic_position(location_hg38)
                if hg38_coords:
                    record["hg38_chromosome"] = hg38_coords.get(
                        "chromosome", chromosome,
                    )
                    if "position" in hg38_coords:
                        record["hg38_position"] = hg38_coords["position"]
                    elif "start" in hg38_coords:
                        record["hg38_start"] = hg38_coords["start"]
                        record["hg38_end"] = hg38_coords["end"]
                    record["assembly_hg38"] = "GRCh38"

            # Parse HG19 coordinates
            if location_hg19:
                hg19_coords = self.parse_genomic_position(location_hg19)
                if hg19_coords:
                    record["hg19_chromosome"] = hg19_coords.get(
                        "chromosome", chromosome,
                    )
                    if "position" in hg19_coords:
                        record["hg19_position"] = hg19_coords["position"]
                    elif "start" in hg19_coords:
                        record["hg19_start"] = hg19_coords["start"]
                        record["hg19_end"] = hg19_coords["end"]
                    record["assembly_hg19"] = "GRCh37"

            # Add gene information
            if gene:
                record["gene"] = gene

            return record

        except Exception as e:
            logger.debug(f"Failed to parse table row: {e}")
            return None
