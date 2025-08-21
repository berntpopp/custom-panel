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

            # Extract rsIDs and coordinates from the list file
            snp_data = self._extract_rsids_from_list(list_path)

            # Create SNP records with coordinates when available
            snps = []
            for data in snp_data:
                rsid = data["rsid"]

                # Extract coordinate information if available
                coord_args = {}
                if "chromosome" in data:
                    coord_args.update(
                        {
                            "chromosome": data["chromosome"],
                            "start": data.get("start"),
                            "end": data.get("end"),
                            "strand": data.get("strand"),
                            "assembly": data.get("assembly"),
                        },
                    )

                snps.append(
                    self.create_snp_record(
                        rsid=rsid,
                        category="identity",
                        panel_specific_name="IDT xGen Human ID Hybridization Panel",
                        **coord_args,
                    ),
                )

            logger.info(
                f"Successfully extracted {len(snps)} SNPs from IDT Hybridization panel",
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

    def _extract_rsids_from_list(self, list_path: Path) -> list[dict[str, Any]]:
        """
        Extract rsIDs with genomic coordinates from the IDT targets list file.

        The .list file is BED-like format with columns:
        chr1	14155402	14155402	+	rs7520386
        - Column 0: Chromosome
        - Column 1: Start position
        - Column 2: End position
        - Column 3: Strand
        - Column 4: rsID

        Args:
            list_path: Path to the list file

        Returns:
            List of dictionaries with rsID and coordinate information
        """
        snp_records = []

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

            logger.info(f"Read {len(df)} rows from {list_path}")

            # Check if we have the expected BED-like format (5 columns)
            if len(df.columns) >= 5:
                logger.info("Found BED-like format with genomic coordinates")

                for _, row in df.iterrows():
                    try:
                        chromosome = str(row.iloc[0]).strip()
                        start = int(row.iloc[1])
                        end = int(row.iloc[2])
                        strand = str(row.iloc[3]).strip()
                        rsid = str(row.iloc[4]).strip()

                        # Validate rsID
                        if self.validate_rsid(rsid):
                            snp_records.append(
                                {
                                    "rsid": rsid,
                                    "chromosome": chromosome,
                                    "start": start,
                                    "end": end,
                                    "strand": strand,
                                    "assembly": "GRCh37",  # IDT typically uses GRCh37/hg19
                                },
                            )
                        else:
                            logger.debug(f"Invalid rsID: {rsid}")

                    except (ValueError, IndexError) as e:
                        logger.debug(f"Skipping malformed row: {e}")
                        continue

            elif len(df.columns) >= 1:
                logger.warning(
                    f"Expected BED format but got {len(df.columns)} columns, extracting rsIDs only",
                )

                # Fallback: try to find rsIDs in any column
                for col in df.columns:
                    col_data = df[col].dropna().astype(str)
                    for value in col_data:
                        if value.startswith("rs") and self.validate_rsid(value):
                            snp_records.append({"rsid": value})

            # If no records found, try line-by-line parsing
            if not snp_records:
                logger.info(
                    "No records found in tabular format, trying line-by-line parsing",
                )

                with open(list_path, encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()

                        # Skip comments and empty lines
                        if line.startswith("@") or not line:
                            continue

                        # Try to parse as tab-delimited BED format
                        parts = line.split("\t")
                        if len(parts) >= 5:
                            try:
                                chromosome = parts[0]
                                start = int(parts[1])
                                end = int(parts[2])
                                strand = parts[3]
                                rsid = parts[4]

                                if self.validate_rsid(rsid):
                                    snp_records.append(
                                        {
                                            "rsid": rsid,
                                            "chromosome": chromosome,
                                            "start": start,
                                            "end": end,
                                            "strand": strand,
                                            "assembly": "GRCh37",
                                        },
                                    )
                            except (ValueError, IndexError):
                                # Extract rsIDs from the line using regex
                                line_rsids = self.extract_rsids_from_text(line)
                                for rsid in line_rsids:
                                    snp_records.append({"rsid": rsid})
                        else:
                            # Extract rsIDs from the line using regex
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
            logger.error(f"Error extracting data from list file: {e}")
            raise
