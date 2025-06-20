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
                        panel_specific_name="Pengelly et al. 2013",
                        **record_kwargs,
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

    def _extract_rsids_from_pdf(self, pdf_path: Path) -> list[dict[str, Any]]:
        """
        Extract rsIDs with genomic coordinates from the Pengelly PDF.

        The PDF contains a table with format:
        Chromosome Position dbSNP rsID Gene Alleles ...

        Args:
            pdf_path: Path to the PDF file

        Returns:
            List of dictionaries with rsID and coordinate information
        """
        snp_records = []

        try:
            with pdfplumber.open(pdf_path) as pdf:
                for page in pdf.pages:
                    # Try table extraction first
                    tables = page.extract_tables()
                    for table in tables:
                        if not table:
                            continue

                        # Look for coordinate table
                        coords_found = self._parse_pengelly_table(table)
                        if coords_found:
                            snp_records.extend(coords_found)
                            logger.info(
                                f"Found Pengelly coordinate table with {len(coords_found)} SNPs"
                            )

                    # If no coordinate table found, fall back to text extraction
                    if not snp_records:
                        text = page.extract_text()
                        if text:
                            # Extract all rsIDs from text
                            rsid_matches = re.findall(
                                r"\b(rs\d+)\b", text, re.IGNORECASE
                            )
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
                f"Extracted {len(unique_records)} SNPs ({with_coords} with coordinates)"
            )
            return unique_records

        except Exception as e:
            logger.error(f"Error extracting rsIDs from PDF: {e}")
            raise

    def _parse_pengelly_table(self, table) -> list[dict[str, Any]]:
        """
        Parse Pengelly table with coordinate information.

        Expected format: Chromosome Position dbSNP rsID Gene Alleles ...

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
            pos_col = None
            rsid_col = None
            gene_col = None
            allele_col = None

            for i, header in enumerate(headers):
                if "chromosome" in header:
                    chr_col = i
                elif "position" in header:
                    pos_col = i
                elif "rsid" in header or "dbsnp" in header:
                    rsid_col = i
                elif "gene" in header:
                    gene_col = i
                elif "allele" in header:
                    allele_col = i

            # If we don't have basic coordinate columns, this isn't the right table
            if chr_col is None or pos_col is None or rsid_col is None:
                return records

            logger.info(
                f"Found Pengelly coordinate table: chr={chr_col}, pos={pos_col}, rsid={rsid_col}"
            )

            # Parse data rows
            for row in table[1:]:  # Skip header
                if not row or len(row) <= max(chr_col, pos_col, rsid_col):
                    continue

                try:
                    chromosome = (
                        str(row[chr_col]).strip()
                        if chr_col < len(row) and row[chr_col]
                        else ""
                    )
                    position_str = (
                        str(row[pos_col]).strip()
                        if pos_col < len(row) and row[pos_col]
                        else ""
                    )
                    rsid = (
                        str(row[rsid_col]).strip()
                        if rsid_col < len(row) and row[rsid_col]
                        else ""
                    )
                    gene = (
                        str(row[gene_col]).strip()
                        if gene_col is not None
                        and gene_col < len(row)
                        and row[gene_col]
                        else ""
                    )
                    alleles = (
                        str(row[allele_col]).strip()
                        if allele_col is not None
                        and allele_col < len(row)
                        and row[allele_col]
                        else ""
                    )

                    # Validate rsID
                    if not self.validate_rsid(rsid):
                        continue

                    record = {"rsid": rsid}

                    # Add chromosome
                    if chromosome:
                        record["chromosome"] = self.clean_chromosome(chromosome)

                    # Add position
                    if position_str and position_str.isdigit():
                        record["position"] = int(position_str)
                        record["assembly"] = "GRCh37"  # Pengelly et al. used GRCh37

                    # Add gene information
                    if gene:
                        record["gene"] = gene

                    # Add allele information
                    if alleles and "/" in alleles:
                        allele_parts = alleles.split("/")
                        if len(allele_parts) == 2:
                            record["ref_allele"] = allele_parts[0].strip()
                            record["alt_allele"] = allele_parts[1].strip()

                    records.append(record)

                except (ValueError, IndexError) as e:
                    logger.debug(f"Could not parse Pengelly table row: {e}")
                    continue

        except Exception as e:
            logger.debug(f"Error parsing Pengelly table: {e}")

        return records
