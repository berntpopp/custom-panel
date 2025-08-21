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
                        panel_specific_name="Illumina AmpliSeq Sample ID Panel",
                        **record_kwargs,
                    ),
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

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[dict[str, Any]]:
        """
        Extract rsIDs with genomic coordinates from the Ampliseq PDF.

        The PDF contains a table with coordinates:
        Code no. Chromosome Start End ID Ref allele Alt allele

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

                    # Try to extract tables first
                    tables = page.extract_tables()
                    for table in tables:
                        if not table:
                            continue

                        # Look for coordinate table
                        coords_found = self._parse_ampliseq_table(table)
                        if coords_found:
                            snp_records.extend(coords_found)
                            logger.info(
                                f"Found coordinate table with {len(coords_found)} SNPs",
                            )

                    # If no coordinate table found, fall back to text extraction
                    if not snp_records:
                        lines = text.split("\n")
                        for line in lines:
                            line = line.strip()

                            # Look for lines starting with number and containing rsID
                            if re.match(r"^\d+.+rs\d+", line):
                                rsid_match = re.search(r"(rs\d+)", line)
                                if rsid_match:
                                    rsid = rsid_match.group(1)
                                    if self.validate_rsid(rsid):
                                        snp_records.append({"rsid": rsid})

                            # Also check for standalone rsIDs
                            elif re.match(r"^rs\d+", line):
                                rsid_match = re.match(r"^(rs\d+)", line)
                                if rsid_match:
                                    rsid = rsid_match.group(1)
                                    if self.validate_rsid(rsid):
                                        snp_records.append({"rsid": rsid})

                            # Extract any rsIDs from the line
                            else:
                                line_rsids = self.extract_rsids_from_text(line)
                                for rsid in line_rsids:
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
                f"Extracted {len(unique_records)} SNPs ({with_coords} with coordinates)",
            )
            return unique_records

        except Exception as e:
            logger.error(f"Error extracting rsIDs from PDF: {e}")
            raise

    def _parse_ampliseq_table(self, table) -> list[dict[str, Any]]:
        """
        Parse AmpliseQ table with coordinate information.

        Expected format: Code no. Chromosome Start End ID Ref allele Alt allele

        Args:
            table: Extracted table from PDF

        Returns:
            List of SNP records with coordinates
        """
        records = []

        try:
            if not table or len(table) < 2:  # Need header + data
                return records

            # Find header row to identify column positions
            header_row = table[0]
            headers = [str(cell).strip().lower() if cell else "" for cell in header_row]

            # Look for coordinate columns
            chr_col = None
            start_col = None
            end_col = None
            id_col = None
            ref_col = None
            alt_col = None

            for i, header in enumerate(headers):
                if "chromosome" in header or header == "chr":
                    chr_col = i
                elif "start" in header:
                    start_col = i
                elif "end" in header:
                    end_col = i
                elif header in ["id", "snp id", "rsid"]:
                    id_col = i
                elif "ref" in header and "allele" in header:
                    ref_col = i
                elif "alt" in header and "allele" in header:
                    alt_col = i

            # If we don't have basic coordinate columns, this isn't the right table
            if chr_col is None or start_col is None or id_col is None:
                return records

            logger.info(
                f"Found AmpliseQ coordinate table: chr={chr_col}, start={start_col}, id={id_col}",
            )

            # Parse data rows
            for row in table[1:]:  # Skip header
                if not row or len(row) <= max(chr_col, start_col, id_col):
                    continue

                try:
                    chromosome = (
                        str(row[chr_col]).strip()
                        if chr_col < len(row) and row[chr_col]
                        else ""
                    )
                    start_str = (
                        str(row[start_col]).strip()
                        if start_col < len(row) and row[start_col]
                        else ""
                    )
                    end_str = (
                        str(row[end_col]).strip()
                        if end_col is not None and end_col < len(row) and row[end_col]
                        else ""
                    )
                    rsid = (
                        str(row[id_col]).strip()
                        if id_col < len(row) and row[id_col]
                        else ""
                    )
                    ref_allele = (
                        str(row[ref_col]).strip()
                        if ref_col is not None and ref_col < len(row) and row[ref_col]
                        else ""
                    )
                    alt_allele = (
                        str(row[alt_col]).strip()
                        if alt_col is not None and alt_col < len(row) and row[alt_col]
                        else ""
                    )

                    # Validate rsID
                    if not self.validate_rsid(rsid):
                        continue

                    record = {"rsid": rsid}

                    # Add chromosome
                    if chromosome:
                        record["chromosome"] = self.clean_chromosome(chromosome)

                    # Add coordinates
                    if start_str and start_str.isdigit():
                        record["start"] = int(start_str)
                        record["assembly"] = "GRCh37"  # Illumina typically uses GRCh37

                    if end_str and end_str.isdigit():
                        record["end"] = int(end_str)

                    # Add allele information
                    if ref_allele:
                        record["ref_allele"] = ref_allele
                    if alt_allele:
                        record["alt_allele"] = alt_allele

                    records.append(record)

                except (ValueError, IndexError) as e:
                    logger.debug(f"Could not parse AmpliseQ table row: {e}")
                    continue

        except Exception as e:
            logger.debug(f"Error parsing AmpliseQ table: {e}")

        return records
