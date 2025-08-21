"""
Identity and ethnicity SNP parsers.

This module contains parsers for identity and ethnicity SNP files
used for sample tracking and ancestry checks.
"""

import json
import logging
import re
from pathlib import Path
from typing import Any

import pandas as pd
import pdfplumber
from bs4 import BeautifulSoup

from .base_snp_parser import BaseSNPParser

logger = logging.getLogger(__name__)


class SimpleRSIDParser(BaseSNPParser):
    """
    Parser for simple text files containing one rsID per line.

    This parser handles plain text files where each line contains
    a single RS identifier, optionally with comments starting with '#'.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a simple text file with one rsID per line.

        Returns:
            DataFrame with columns: rsid, source, category
        """
        self.validate_file_exists()

        logger.info(f"Parsing simple rsID file: {self.file_path}")

        try:
            # Read file and filter out comments and empty lines
            with open(self.file_path, encoding="utf-8") as f:
                lines = f.readlines()

            rsids = []
            for line in lines:
                line = line.strip()
                # Skip empty lines and comments
                if line and not line.startswith("#"):
                    rsids.append(line)

            if not rsids:
                logger.warning(f"No valid rsIDs found in {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Create DataFrame
            df = pd.DataFrame(
                {"rsid": rsids, "source": self.name, "category": "identity"},
            )

            logger.info(f"Successfully parsed {len(df)} rsIDs from {self.file_path}")

            # Validate and standardize output
            return self.validate_output_format(df)

        except Exception as e:
            logger.error(f"Error parsing simple rsID file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse simple rsID file: {e}") from e


class ExcelRSIDParser(BaseSNPParser):
    """
    Parser for Excel files containing rsIDs and source information.

    This parser handles .xlsx files with configurable column mappings
    for rsID and source group information.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse an Excel file with rsID and source columns.

        Returns:
            DataFrame with columns: rsid, source, category
        """
        self.validate_file_exists()

        logger.info(f"Parsing Excel rsID file: {self.file_path}")

        # Get column mappings from config
        rsid_column = self.config.get("rsid_column", "rs_id")
        source_column = self.config.get("source_column", "group")
        sheet_name = self.config.get("sheet_name", 0)  # Default to first sheet

        try:
            # Read Excel file
            df = pd.read_excel(self.file_path, sheet_name=sheet_name)

            # Check for required columns
            if rsid_column not in df.columns:
                raise ValueError(
                    f"Required rsID column '{rsid_column}' not found. Available columns: {list(df.columns)}",
                )

            # Extract rsIDs
            rsids = df[rsid_column].dropna().tolist()

            if not rsids:
                logger.warning(
                    f"No valid rsIDs found in column '{rsid_column}' of {self.file_path}",
                )
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Create base DataFrame
            result_df = pd.DataFrame(
                {
                    "rsid": rsids,
                    "source": self.name,
                    "category": "ethnicity",  # Excel files are typically for ethnicity
                },
            )

            # Add source group information if available
            if source_column in df.columns:
                source_groups = df[source_column].dropna().tolist()
                if len(source_groups) == len(rsids):
                    result_df["source_group"] = source_groups
                else:
                    logger.warning(
                        f"Source column '{source_column}' length doesn't match rsID column",
                    )

            logger.info(
                f"Successfully parsed {len(result_df)} rsIDs from {self.file_path}",
            )

            # Validate and standardize output
            return self.validate_output_format(result_df)

        except Exception as e:
            logger.error(f"Error parsing Excel rsID file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse Excel rsID file: {e}") from e


class PDFRSIDParser(BaseSNPParser):
    """
    Parser for PDF files containing rsIDs.

    This parser extracts text from PDF files and searches for rsID patterns.
    Supports different extraction patterns based on panel type.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a PDF file for rsIDs using configurable patterns.

        Returns:
            DataFrame with columns: rsid, source, category
        """
        self.validate_file_exists()

        logger.info(f"Parsing PDF rsID file: {self.file_path}")

        # Get pattern configuration
        pattern_type = self.config.get("pattern_type", "default")

        try:
            rsids = []

            with pdfplumber.open(self.file_path) as pdf:
                for page in pdf.pages:
                    text = page.extract_text()
                    if not text:
                        continue

                    # Split into lines and process based on pattern type
                    lines = text.split("\n")

                    if pattern_type == "pengelly":
                        # Pattern: ^rs[0-9]* at start of line
                        for line in lines:
                            line = line.strip()
                            if re.match(r"^rs\d+", line):
                                rsid_match = re.match(r"^(rs\d+)", line)
                                if rsid_match:
                                    rsids.append(rsid_match.group(1))

                    elif pattern_type == "eurogentest":
                        # Pattern: chr.+rs[0-9]* (chromosome info followed by rsID)
                        for line in lines:
                            if re.search(r"chr.+rs\d+", line):
                                # Extract rsID from line with chr info
                                parts = line.split()
                                for part in parts:
                                    if re.match(r"^rs\d+$", part):
                                        rsids.append(part)

                    elif pattern_type == "ampliseq":
                        # Pattern: ^[0-9].+rs[0-9]* (number at start, rsID somewhere)
                        for line in lines:
                            if re.match(r"^\d+.+rs\d+", line):
                                # Extract rsID from line
                                rsid_match = re.search(r"(rs\d+)", line)
                                if rsid_match:
                                    rsids.append(rsid_match.group(1))

                    else:  # default pattern
                        # Generic rsID extraction
                        rsid_matches = re.findall(r"\b(rs\d+)\b", text)
                        rsids.extend(rsid_matches)

            # Remove duplicates while preserving order
            seen = set()
            unique_rsids = []
            for rsid in rsids:
                if rsid not in seen:
                    seen.add(rsid)
                    unique_rsids.append(rsid)

            if not unique_rsids:
                logger.warning(f"No valid rsIDs found in PDF {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Create DataFrame
            df = pd.DataFrame(
                {"rsid": unique_rsids, "source": self.name, "category": "identity"},
            )

            logger.info(
                f"Successfully parsed {len(df)} unique rsIDs from PDF {self.file_path}",
            )

            # Validate and standardize output
            return self.validate_output_format(df)

        except Exception as e:
            logger.error(f"Error parsing PDF file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse PDF file: {e}") from e


class HTMLRSIDParser(BaseSNPParser):
    """
    Parser for HTML files containing rsIDs.

    This parser extracts rsIDs from HTML tables or other structured elements.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse an HTML file for rsIDs in tables.

        Returns:
            DataFrame with columns: rsid, source, category
        """
        self.validate_file_exists()

        logger.info(f"Parsing HTML rsID file: {self.file_path}")

        # Get configuration
        table_class = self.config.get("table_class", "table-condensed")
        rsid_column = self.config.get("rsid_column", "SNP ID")

        try:
            with open(self.file_path, encoding="utf-8") as f:
                html_content = f.read()

            soup = BeautifulSoup(html_content, "html.parser")

            rsids = []

            # Find tables with specified class
            tables = soup.find_all("table", class_=table_class)

            if not tables:
                # Try finding any table
                tables = soup.find_all("table")

            for table in tables:
                # Try to parse as pandas DataFrame
                try:
                    df_list = pd.read_html(str(table))
                    if df_list:
                        df = df_list[0]

                        # Look for rsID column
                        if rsid_column in df.columns:
                            table_rsids = df[rsid_column].dropna().tolist()
                            rsids.extend(
                                [
                                    str(rsid)
                                    for rsid in table_rsids
                                    if str(rsid).startswith("rs")
                                ],
                            )
                        else:
                            # Search all columns for rsIDs
                            for col in df.columns:
                                col_data = df[col].dropna().astype(str)
                                col_rsids = col_data[
                                    col_data.str.match(r"^rs\d+$", na=False)
                                ]
                                rsids.extend(col_rsids.tolist())
                except Exception as e:
                    logger.warning(f"Failed to parse table: {e}")

            # If no tables found or no rsIDs in tables, search entire HTML
            if not rsids:
                rsid_matches = re.findall(r"\b(rs\d+)\b", html_content)
                rsids = list(set(rsid_matches))  # Remove duplicates

            if not rsids:
                logger.warning(f"No valid rsIDs found in HTML {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Create DataFrame
            df = pd.DataFrame(
                {"rsid": rsids, "source": self.name, "category": "identity"},
            )

            logger.info(
                f"Successfully parsed {len(df)} rsIDs from HTML {self.file_path}",
            )

            # Validate and standardize output
            return self.validate_output_format(df)

        except Exception as e:
            logger.error(f"Error parsing HTML file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse HTML file: {e}") from e


class TSVListParser(BaseSNPParser):
    """
    Parser for TSV/List files containing rsIDs.

    This parser handles tab-separated or list files with configurable
    column positions for rsID extraction.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a TSV/List file for rsIDs.

        Returns:
            DataFrame with columns: rsid, source, category
        """
        self.validate_file_exists()

        logger.info(f"Parsing TSV/List rsID file: {self.file_path}")

        # Get configuration
        rsid_column = self.config.get("rsid_column", 4)  # 0-based index or name
        comment_char = self.config.get("comment", "@")
        delimiter = self.config.get("delimiter", "\t")

        try:
            # Read file with pandas
            if isinstance(rsid_column, int):
                # Column index provided
                df = pd.read_csv(
                    self.file_path,
                    sep=delimiter,
                    comment=comment_char,
                    header=None,
                    engine="python",
                )

                if rsid_column >= len(df.columns):
                    raise ValueError(
                        f"Column index {rsid_column} out of range. File has {len(df.columns)} columns.",
                    )

                rsid_data = df.iloc[:, rsid_column]
            else:
                # Column name provided
                df = pd.read_csv(
                    self.file_path, sep=delimiter, comment=comment_char, engine="python",
                )

                if rsid_column not in df.columns:
                    raise ValueError(
                        f"Column '{rsid_column}' not found. Available: {list(df.columns)}",
                    )

                rsid_data = df[rsid_column]

            # Extract rsIDs
            rsids = []
            for value in rsid_data.dropna():
                value_str = str(value).strip()
                if re.match(r"^rs\d+$", value_str):
                    rsids.append(value_str)

            if not rsids:
                logger.warning(f"No valid rsIDs found in {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Create DataFrame
            result_df = pd.DataFrame(
                {"rsid": rsids, "source": self.name, "category": "identity"},
            )

            logger.info(
                f"Successfully parsed {len(result_df)} rsIDs from {self.file_path}",
            )

            # Validate and standardize output
            return self.validate_output_format(result_df)

        except Exception as e:
            logger.error(f"Error parsing TSV/List file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse TSV/List file: {e}") from e


class ScrapedJSONParser(BaseSNPParser):
    """
    Parser for scraped JSON files from SNP scrapers.

    This parser handles JSON files that contain SNP data in the format
    output by the custom SNP scrapers (pengelly, eurogentest, etc.).
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a scraped JSON file for rsIDs and metadata.

        Returns:
            DataFrame with columns: rsid, source, category, plus any additional metadata
        """
        self.validate_file_exists()

        logger.info(f"Parsing scraped JSON file: {self.file_path}")

        try:
            # Read JSON file
            with open(self.file_path, encoding="utf-8") as f:
                data = json.load(f)

            # Extract SNPs from the JSON structure
            snps = data.get("snps", [])

            if not snps:
                logger.warning(f"No SNPs found in JSON file {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])

            # Convert to DataFrame
            df = pd.DataFrame(snps)

            # Ensure required columns exist
            if "rsid" not in df.columns:
                raise ValueError(
                    f"Required 'rsid' column not found in {self.file_path}",
                )

            # Rename 'snp' column to 'rsid' if it exists
            if "snp" in df.columns and "rsid" not in df.columns:
                df = df.rename(columns={"snp": "rsid"})

            # Set default values for missing columns
            if "source" not in df.columns:
                df["source"] = self.name
            if "category" not in df.columns:
                df["category"] = "identity"

            # Keep only non-null rsIDs
            df = df.dropna(subset=["rsid"])

            if df.empty:
                logger.warning(f"No valid rsIDs found in JSON file {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])

            logger.info(
                f"Successfully parsed {len(df)} SNPs from scraped JSON {self.file_path}",
            )

            # Validate and standardize output
            return self.validate_output_format(df)

        except Exception as e:
            logger.error(f"Error parsing scraped JSON file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse scraped JSON file: {e}") from e


def create_identity_parser(file_path: Path, config: dict[str, Any]) -> BaseSNPParser:
    """
    Factory function to create the appropriate identity parser based on configuration.

    Args:
        file_path: Path to the SNP file
        config: Parser configuration

    Returns:
        Appropriate parser instance

    Raises:
        ValueError: If parser type is not supported
    """
    parser_type = config.get("parser", "simple_rsid")

    if parser_type == "simple_rsid":
        return SimpleRSIDParser(file_path, config)
    elif parser_type == "excel_rsid":
        return ExcelRSIDParser(file_path, config)
    elif parser_type == "pdf_rsid":
        return PDFRSIDParser(file_path, config)
    elif parser_type == "html_rsid":
        return HTMLRSIDParser(file_path, config)
    elif parser_type == "tsv_list":
        return TSVListParser(file_path, config)
    elif parser_type == "scraped_json":
        return ScrapedJSONParser(file_path, config)
    else:
        raise ValueError(f"Unsupported identity parser type: {parser_type}")
