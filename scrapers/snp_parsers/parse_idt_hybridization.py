"""
IDT xGen Human ID Hybridization panel SNP scraper.

Extracts identity SNPs from the IDT targets list file for their hybridization-based panel.
Reference: https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/supplementary-product-info/xgen-human-id-hyb-panel-targetsfbea8d1532796e2eaa53ff00001c1b3c.list?sfvrsn=3693e307_2&download=true
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from custom_panel.sources_snp.downloader import PanelDownloader

from .base_snp_scraper import BaseSNPScraper

logger = logging.getLogger(__name__)


class IDTHybridizationParser(BaseSNPScraper):
    """
    Parser for IDT xGen Human ID Hybridization panel.

    This parser extracts rsIDs from the IDT targets list file (.list format)
    which contains target information for their hybridization capture panel.
    """

    def parse(self) -> dict[str, Any]:
        """
        Parse the IDT targets list file and extract SNP data.

        Returns:
            Dictionary containing panel data and SNPs
        """
        logger.info(f"Parsing IDT Hybridization panel from {self.url}")

        try:
            # Download list file
            with PanelDownloader() as downloader:
                list_path = downloader.download(self.url, "list", "requests")

            # Extract rsIDs from the list file
            rsids = self._extract_rsids_from_list(list_path)

            # Create SNP records
            snps = []
            for rsid in rsids:
                snps.append(
                    self.create_snp_record(
                        rsid=rsid,
                        category="identity",
                        panel_specific_name="IDT xGen Human ID Hybridization Panel",
                    )
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from IDT Hybridization panel"
            )

            return {
                "panel_name": "idt_xgen_human_id_hybridization_panel",
                "source_url": self.url,
                "description": "IDT xGen Human ID Hybridization Panel - Capture-based sample tracking",
                "snps": snps,
                "metadata": {
                    "vendor": "Integrated DNA Technologies",
                    "technology": "Hybridization capture",
                    "source_url": self.url,
                    "snp_count": len(snps),
                },
            }

        except Exception as e:
            logger.error(f"Error parsing IDT Hybridization panel: {e}")
            raise

    def _extract_rsids_from_list(self, list_path: Path) -> list[str]:
        """
        Extract rsIDs from the IDT targets list file.

        The .list file is typically tab-delimited with multiple columns.
        Based on the R script, rsIDs are in column 5 (0-based index 4).
        Lines starting with @ are comments.

        Args:
            list_path: Path to the list file

        Returns:
            List of unique rsIDs
        """
        rsids = []

        try:
            # Read the file as tab-delimited
            df = pd.read_csv(
                list_path,
                sep="\t",
                comment="@",
                header=None,
                engine="python",
                on_bad_lines="skip",  # Skip malformed lines
            )

            # Check if we have at least 5 columns (0-based index 4)
            if len(df.columns) < 5:
                logger.warning(f"Expected at least 5 columns, got {len(df.columns)}")
                # Try to find rsIDs in any column
                for col in df.columns:
                    col_data = df[col].dropna().astype(str)
                    for value in col_data:
                        if value.startswith("rs") and self.validate_rsid(value):
                            rsids.append(value)
            else:
                # Extract from column 5 (index 4) as per R script
                rsid_column = df.iloc[:, 4]  # 5th column (0-based index 4)

                for value in rsid_column.dropna():
                    value_str = str(value).strip()
                    if value_str.startswith("rs") and self.validate_rsid(value_str):
                        rsids.append(value_str)

            # Also try parsing as a more flexible format
            if not rsids:
                logger.info(
                    "No rsIDs found in expected column, trying alternative parsing"
                )

                # Read line by line and extract rsIDs
                with open(list_path, encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()

                        # Skip comments
                        if line.startswith("@") or not line:
                            continue

                        # Extract rsIDs from the line
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
            logger.error(f"Error extracting rsIDs from list file: {e}")
            raise
