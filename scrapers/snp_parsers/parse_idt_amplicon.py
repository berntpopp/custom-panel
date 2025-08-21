"""
IDT xGen Sample ID Amplicon panel SNP scraper.

Extracts identity SNPs from the IDT website for their amplicon-based sample ID panel.
Reference: https://www.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels/sample-id-amp-panel
"""

import logging
from typing import Any

from bs4 import BeautifulSoup

from .base_snp_scraper import BaseSNPScraper

logger = logging.getLogger(__name__)


class IDTAmpliconParser(BaseSNPScraper):
    """
    Parser for IDT xGen Sample ID Amplicon panel.

    This parser extracts rsIDs from the IDT website's HTML page
    that lists SNPs in their sample identification amplicon panel.
    """

    def parse(self) -> dict[str, Any]:
        """
        Parse the IDT website and extract SNP data.

        Returns:
            Dictionary containing panel data and SNPs
        """
        logger.info(f"Parsing IDT Amplicon panel from {self.url}")

        try:
            # Try Selenium first for JavaScript rendering
            html_content = None
            try:
                html_content = self.fetch_content(use_selenium=True)
                logger.info("Successfully fetched content with Selenium")
            except Exception as selenium_error:
                logger.warning(f"Selenium failed: {selenium_error}")
                logger.info("Falling back to requests method")
                try:
                    html_content = self.fetch_content(use_selenium=False)
                    logger.info("Successfully fetched content with requests")
                except Exception as requests_error:
                    logger.error(
                        f"Both Selenium and requests failed. Selenium: {selenium_error}, Requests: {requests_error}",
                    )
                    raise Exception(
                        f"All download methods failed. Last error: {requests_error}",
                    ) from requests_error

            if not html_content:
                raise Exception("No content retrieved from any download method")

            # Parse HTML
            soup = BeautifulSoup(html_content, "html.parser")

            # Extract rsIDs and coordinates from the page
            snp_data = self._extract_rsids_from_html(soup)

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
                        panel_specific_name="IDT xGen Sample ID Amplicon Panel",
                        **record_kwargs,
                    ),
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from IDT Amplicon panel",
            )

            return {
                "panel_name": "idt_xgen_sample_id_amplicon_panel",
                "source_url": self.url,
                "description": "IDT xGen Sample ID Amplicon Panel - PCR-based sample tracking",
                "snps": snps,
                "metadata": {
                    "vendor": "Integrated DNA Technologies",
                    "technology": "Amplicon sequencing",
                    "source_url": self.url,
                    "snp_count": len(snps),
                },
            }

        except Exception as e:
            logger.error(f"Error parsing IDT Amplicon panel: {e}")
            raise

    def _extract_rsids_from_html(self, soup: BeautifulSoup) -> list[dict[str, Any]]:
        """
        Extract rsIDs with genomic coordinates from the IDT HTML page.

        The page contains a table with columns:
        Chr    POS    SNP ID    REF    ALT    [sample columns...]

        Args:
            soup: BeautifulSoup object of the page

        Returns:
            List of dictionaries with rsID and coordinate information
        """
        snp_records = []

        # Try multiple strategies to find the SNP data

        # Strategy 1: Look for tables with coordinate columns
        tables = soup.find_all("table")
        for table in tables:
            records = self._parse_coordinate_table(table)
            if records:
                snp_records.extend(records)
                logger.info(f"Found coordinate table with {len(records)} SNPs")
                break

        # Strategy 2: Look for tables with specific classes
        if not snp_records:
            table_classes = ["table-condensed", "table", "data-table", "snp-table"]
            rsid_columns = ["SNP ID", "rsID", "RS ID", "SNP", "Marker"]

            for table_class in table_classes:
                tables = soup.find_all("table", class_=table_class)
                if tables:
                    for rsid_col in rsid_columns:
                        table_rsids = self.parse_table_for_rsids(
                            soup,
                            f"table.{table_class}",
                            rsid_col,
                        )
                        if table_rsids:
                            snp_records.extend([{"rsid": rsid} for rsid in table_rsids])
                            break
                    if snp_records:
                        break

        # Strategy 3: Look for any table and check all columns
        if not snp_records:
            table_rsids = self.parse_table_for_rsids(soup)
            if table_rsids:
                snp_records.extend([{"rsid": rsid} for rsid in table_rsids])

        # Strategy 4: Look for rsIDs in lists or divs
        if not snp_records:
            # Check for rsIDs in list items
            for li in soup.find_all("li"):
                text = li.get_text()
                rsids = self.extract_rsids_from_text(text)
                snp_records.extend([{"rsid": rsid} for rsid in rsids])

            # Check for rsIDs in divs with specific classes
            for div_class in ["snp-list", "marker-list", "panel-content"]:
                divs = soup.find_all("div", class_=div_class)
                for div in divs:
                    text = div.get_text()
                    rsids = self.extract_rsids_from_text(text)
                    snp_records.extend([{"rsid": rsid} for rsid in rsids])

        # Strategy 5: Extract from entire page text as last resort
        if not snp_records:
            page_text = soup.get_text()
            rsids = self.extract_rsids_from_text(page_text)
            snp_records.extend([{"rsid": rsid} for rsid in rsids])

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

    def _parse_coordinate_table(self, table) -> list[dict[str, Any]]:
        """
        Parse a table with coordinate information.

        Expected format: Chr    POS    SNP ID    REF    ALT    [sample columns...]

        Args:
            table: BeautifulSoup table element

        Returns:
            List of SNP records with coordinates
        """
        records = []

        try:
            rows = table.find_all("tr")
            if not rows:
                return records

            # Find header row to identify column positions
            header_row = rows[0]
            headers = [
                th.get_text().strip() for th in header_row.find_all(["th", "td"])
            ]

            # Look for coordinate columns
            chr_col = None
            pos_col = None
            snp_col = None
            ref_col = None
            alt_col = None

            for i, header in enumerate(headers):
                header_lower = header.lower()
                if header_lower in ["chr", "chromosome"]:
                    chr_col = i
                elif header_lower in ["pos", "position"]:
                    pos_col = i
                elif header_lower in ["snp id", "snp_id", "rsid", "snp"]:
                    snp_col = i
                elif header_lower in ["ref", "reference"]:
                    ref_col = i
                elif header_lower in ["alt", "alternative"]:
                    alt_col = i

            # If we don't have basic coordinate columns, this isn't the right table
            if chr_col is None or pos_col is None or snp_col is None:
                return records

            logger.info(
                f"Found coordinate table: chr={chr_col}, pos={pos_col}, snp={snp_col}",
            )

            # Parse data rows
            for row in rows[1:]:  # Skip header
                cells = [td.get_text().strip() for td in row.find_all(["td", "th"])]
                if len(cells) <= max(chr_col, pos_col, snp_col):
                    continue

                try:
                    chromosome = cells[chr_col] if chr_col < len(cells) else ""
                    position_str = cells[pos_col] if pos_col < len(cells) else ""
                    rsid = cells[snp_col] if snp_col < len(cells) else ""
                    ref_allele = (
                        cells[ref_col]
                        if ref_col is not None and ref_col < len(cells)
                        else ""
                    )
                    alt_allele = (
                        cells[alt_col]
                        if alt_col is not None and alt_col < len(cells)
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
                        record["assembly"] = "GRCh37"  # IDT typically uses GRCh37

                    # Add allele information
                    if ref_allele:
                        record["ref_allele"] = ref_allele
                    if alt_allele:
                        record["alt_allele"] = alt_allele

                    records.append(record)

                except (ValueError, IndexError) as e:
                    logger.debug(f"Could not parse table row: {e}")
                    continue

        except Exception as e:
            logger.debug(f"Error parsing coordinate table: {e}")

        return records
