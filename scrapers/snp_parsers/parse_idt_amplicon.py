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
                        f"Both Selenium and requests failed. Selenium: {selenium_error}, Requests: {requests_error}"
                    )
                    raise Exception(
                        f"All download methods failed. Last error: {requests_error}"
                    ) from requests_error

            if not html_content:
                raise Exception("No content retrieved from any download method")

            # Parse HTML
            soup = BeautifulSoup(html_content, "html.parser")

            # Extract rsIDs from the page
            rsids = self._extract_rsids_from_html(soup)

            # Create SNP records
            snps = []
            for rsid in rsids:
                snps.append(
                    self.create_snp_record(
                        rsid=rsid,
                        category="identity",
                        panel_specific_name="IDT xGen Sample ID Amplicon Panel",
                    )
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from IDT Amplicon panel"
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

    def _extract_rsids_from_html(self, soup: BeautifulSoup) -> list[str]:
        """
        Extract rsIDs from the IDT HTML page.

        The page typically contains a table with class "table-condensed" or similar
        with a column "SNP ID" containing the rsIDs.

        Args:
            soup: BeautifulSoup object of the page

        Returns:
            List of unique rsIDs
        """
        rsids = []

        # Try multiple strategies to find the SNP data

        # Strategy 1: Look for tables with specific classes
        table_classes = ["table-condensed", "table", "data-table", "snp-table"]
        rsid_columns = ["SNP ID", "rsID", "RS ID", "SNP", "Marker"]

        for table_class in table_classes:
            tables = soup.find_all("table", class_=table_class)
            if tables:
                for rsid_col in rsid_columns:
                    table_rsids = self.parse_table_for_rsids(
                        soup, f"table.{table_class}", rsid_col
                    )
                    if table_rsids:
                        rsids.extend(table_rsids)
                        break
                if rsids:
                    break

        # Strategy 2: Look for any table and check all columns
        if not rsids:
            rsids = self.parse_table_for_rsids(soup)

        # Strategy 3: Look for rsIDs in lists or divs
        if not rsids:
            # Check for rsIDs in list items
            for li in soup.find_all("li"):
                text = li.get_text()
                rsids.extend(self.extract_rsids_from_text(text))

            # Check for rsIDs in divs with specific classes
            for div_class in ["snp-list", "marker-list", "panel-content"]:
                divs = soup.find_all("div", class_=div_class)
                for div in divs:
                    text = div.get_text()
                    rsids.extend(self.extract_rsids_from_text(text))

        # Strategy 4: Extract from entire page text as last resort
        if not rsids:
            page_text = soup.get_text()
            rsids = self.extract_rsids_from_text(page_text)

        # Remove duplicates while preserving order
        seen = set()
        unique_rsids = []
        for rsid in rsids:
            if rsid not in seen:
                seen.add(rsid)
                unique_rsids.append(rsid)

        return unique_rsids
