"""
Polygenic Risk Score (PRS) SNP parsers.

This module contains parsers for PRS SNP files which typically have
more complex formats with additional metadata like effect alleles,
weights, and chromosomal positions.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

from .base_snp_parser import BaseSNPParser

logger = logging.getLogger(__name__)


class BCACParser(BaseSNPParser):
    """
    Parser for BCAC (Breast Cancer Association Consortium) PRS files.
    
    This parser handles the specific format used by BCAC PRS files,
    which typically contain columns for rsID, chromosome, position,
    effect allele, and other PRS-specific metadata.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a BCAC PRS CSV file.
        
        Returns:
            DataFrame with columns: rsid, source, category, plus PRS metadata
        """
        self.validate_file_exists()
        
        logger.info(f"Parsing BCAC PRS file: {self.file_path}")
        
        # Get column mappings from config
        rsid_column = self.config.get("rsid_column", "hm_rsid")
        chr_column = self.config.get("chr_column", "hm_chr")
        pos_column = self.config.get("pos_column", "hm_pos")
        effect_allele_column = self.config.get("effect_allele_column", "effect_allele")
        
        try:
            # Read CSV file
            df = pd.read_csv(self.file_path)
            
            # Check for required rsID column
            if rsid_column not in df.columns:
                raise ValueError(f"Required rsID column '{rsid_column}' not found. Available columns: {list(df.columns)}")
            
            # Extract rsIDs and remove empty ones
            rsids = df[rsid_column].dropna()
            if rsids.empty:
                logger.warning(f"No valid rsIDs found in column '{rsid_column}' of {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])
            
            # Create base DataFrame with required columns
            result_df = pd.DataFrame({
                "rsid": rsids,
                "source": self.name,
                "category": "prs"
            })
            
            # Add PRS-specific metadata columns if available
            metadata_columns = {
                "chromosome": chr_column,
                "position": pos_column,
                "effect_allele": effect_allele_column,
            }
            
            for result_col, source_col in metadata_columns.items():
                if source_col in df.columns:
                    # Align with rsid index to ensure proper mapping
                    result_df[result_col] = df.loc[rsids.index, source_col].values
                else:
                    logger.debug(f"Optional column '{source_col}' not found in {self.file_path}")
            
            # Add any additional columns that might be useful for PRS
            additional_columns = ["effect_weight", "OR", "beta", "se", "pvalue", "freq"]
            for col in additional_columns:
                if col in df.columns:
                    result_df[col] = df.loc[rsids.index, col].values
            
            logger.info(f"Successfully parsed {len(result_df)} PRS SNPs from {self.file_path}")
            
            # Validate and standardize output
            return self.validate_output_format(result_df)
            
        except Exception as e:
            logger.error(f"Error parsing BCAC PRS file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse BCAC PRS file: {e}") from e


class GenericPRSParser(BaseSNPParser):
    """
    Generic parser for PRS files with configurable column mappings.
    
    This parser can handle various PRS file formats by allowing
    flexible column mapping configuration.
    """

    def parse(self) -> pd.DataFrame:
        """
        Parse a generic PRS file with configurable column mappings.
        
        Returns:
            DataFrame with columns: rsid, source, category, plus metadata
        """
        self.validate_file_exists()
        
        logger.info(f"Parsing generic PRS file: {self.file_path}")
        
        # Get required column mapping
        rsid_column = self.config.get("rsid_column")
        if not rsid_column:
            raise ValueError("rsid_column must be specified in config for generic PRS parser")
        
        # Get file format (default to CSV)
        file_format = self.config.get("format", "csv")
        separator = self.config.get("separator", ",")
        
        try:
            # Read file based on format
            if file_format.lower() == "csv":
                df = pd.read_csv(self.file_path, sep=separator)
            elif file_format.lower() in ["xlsx", "excel"]:
                sheet_name = self.config.get("sheet_name", 0)
                df = pd.read_excel(self.file_path, sheet_name=sheet_name)
            elif file_format.lower() == "tsv":
                df = pd.read_csv(self.file_path, sep="\t")
            else:
                raise ValueError(f"Unsupported file format: {file_format}")
            
            # Check for required rsID column
            if rsid_column not in df.columns:
                raise ValueError(f"Required rsID column '{rsid_column}' not found. Available columns: {list(df.columns)}")
            
            # Extract rsIDs
            rsids = df[rsid_column].dropna()
            if rsids.empty:
                logger.warning(f"No valid rsIDs found in column '{rsid_column}' of {self.file_path}")
                return pd.DataFrame(columns=["rsid", "source", "category"])
            
            # Create base DataFrame
            result_df = pd.DataFrame({
                "rsid": rsids,
                "source": self.name,
                "category": "prs"
            })
            
            # Add any additional columns specified in column mappings
            column_mappings = self.config.get("column_mappings", {})
            for result_col, source_col in column_mappings.items():
                if source_col in df.columns:
                    result_df[result_col] = df.loc[rsids.index, source_col].values
            
            logger.info(f"Successfully parsed {len(result_df)} PRS SNPs from {self.file_path}")
            
            # Validate and standardize output
            return self.validate_output_format(result_df)
            
        except Exception as e:
            logger.error(f"Error parsing generic PRS file {self.file_path}: {e}")
            raise ValueError(f"Failed to parse generic PRS file: {e}") from e


def create_prs_parser(file_path: Path, config: dict[str, Any]) -> BaseSNPParser:
    """
    Factory function to create the appropriate PRS parser based on configuration.
    
    Args:
        file_path: Path to the SNP file
        config: Parser configuration
        
    Returns:
        Appropriate parser instance
        
    Raises:
        ValueError: If parser type is not supported
    """
    parser_type = config.get("parser", "bcac_prs")
    
    if parser_type == "bcac_prs":
        return BCACParser(file_path, config)
    elif parser_type == "generic_prs":
        return GenericPRSParser(file_path, config)
    else:
        raise ValueError(f"Unsupported PRS parser type: {parser_type}")