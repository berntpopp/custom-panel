"""
Somatic commercial panels data source extractor.

This module processes somatic cancer panels from commercial providers that offer
direct file downloads (PDF, Excel, ZIP formats) as identified in the original R implementation.
"""

import logging
import re
import tempfile
import zipfile
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pandas as pd
import requests
from pypdf import PdfReader

from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)


def fetch_somatic_commercial_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch somatic commercial panel data from configured sources.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with somatic commercial panel data
    """
    somatic_config = config.get("data_sources", {}).get("somatic_commercial", {})

    if not somatic_config.get("enabled", True):
        logger.info("Somatic commercial panels data source is disabled")
        return pd.DataFrame()

    panels_config = somatic_config.get("panels", [])
    if not panels_config:
        logger.warning("No somatic commercial panels configured")
        return pd.DataFrame()

    all_dataframes = []

    for panel_config in panels_config:
        panel_name = panel_config.get("name", "Unknown_Panel")
        logger.info(f"Processing somatic commercial panel: {panel_name}")

        try:
            df = process_somatic_panel(panel_config)
            if df is not None and not df.empty:
                all_dataframes.append(df)
                logger.info(f"Successfully processed {len(df)} genes from {panel_name}")
            else:
                logger.warning(f"No genes extracted from {panel_name}")
        except Exception as e:
            logger.error(f"Error processing panel {panel_name}: {e}")
            continue

    # Combine all dataframes
    if all_dataframes:
        combined_df = pd.concat(all_dataframes, ignore_index=True)
        logger.info(
            f"Fetched {len(combined_df)} total gene records from somatic commercial panels"
        )
        return combined_df
    else:
        logger.warning("No somatic commercial panel data was successfully processed")
        return pd.DataFrame()


def process_somatic_panel(panel_config: dict[str, Any]) -> pd.DataFrame | None:
    """
    Process a single somatic commercial panel.

    Args:
        panel_config: Panel configuration dictionary

    Returns:
        Standardized DataFrame or None if processing failed
    """
    panel_name = panel_config.get("name", "Unknown_Panel")
    url = panel_config.get("url")
    evidence_score = panel_config.get("evidence_score", 0.8)
    category = panel_config.get("category", "somatic")  # Default to somatic
    file_type = panel_config.get("type", "auto")  # pdf, xlsx, zip, auto

    if not url:
        logger.error(f"No URL specified for panel: {panel_name}")
        return None

    # Download file
    try:
        file_path = download_panel_file(url, panel_name, file_type)
        if not file_path:
            return None
    except Exception as e:
        logger.error(f"Error downloading file for {panel_name}: {e}")
        return None

    # Process based on panel name (matching R implementation logic)
    try:
        if "idtxgen_pan_cancer" in panel_name.lower():
            genes = process_idtxgen_panel(file_path)
        elif "illumina_trusight_oncology" in panel_name.lower():
            genes = process_illumina_trusight_panel(file_path)
        elif "sureselect_cancer_all_in_one" in panel_name.lower():
            genes = process_sureselect_allinone_panel(file_path)
        elif "sureselect_genomics_glasgow" in panel_name.lower():
            genes = process_sureselect_glasgow_panel(file_path, panel_name)
        elif "ion_ampliseq_cancer" in panel_name.lower():
            genes = process_ion_ampliseq_panel(file_path)
        elif "illumina_comprehensive_panel" in panel_name.lower():
            genes = process_illumina_comprehensive_panel(file_path)
        else:
            logger.warning(
                f"Unknown panel type for {panel_name}, using generic processing"
            )
            genes = process_generic_panel(file_path, file_type)

        if not genes:
            logger.warning(f"No genes found in {panel_name}")
            return None

        # Create standardized dataframe
        evidence_scores = [evidence_score] * len(genes)
        source_details = [f"URL:{url}|Category:{category}"] * len(genes)

        df = create_standard_dataframe(
            genes=genes,
            source_name=panel_name,
            evidence_scores=evidence_scores,
            source_details=source_details,
            gene_names_reported=genes,
        )

        # Add temporary category column for use by the merger
        df["category"] = category

        return df

    except Exception as e:
        logger.error(f"Error processing {panel_name}: {e}")
        return None
    finally:
        # Clean up downloaded file
        if file_path and file_path.exists():
            file_path.unlink()


def download_panel_file(url: str, panel_name: str, file_type: str) -> Path | None:
    """
    Download a panel file from URL.

    Args:
        url: URL to download from
        panel_name: Name of the panel for logging
        file_type: Expected file type (pdf, xlsx, zip, auto)

    Returns:
        Path to downloaded file or None if failed
    """
    try:
        # Determine file extension
        if file_type == "auto":
            parsed_url = urlparse(url)
            file_ext = Path(parsed_url.path).suffix.lower()
            if not file_ext:
                # Try to guess from URL
                if "pdf" in url.lower():
                    file_ext = ".pdf"
                elif "xlsx" in url.lower() or "excel" in url.lower():
                    file_ext = ".xlsx"
                elif "zip" in url.lower():
                    file_ext = ".zip"
                else:
                    file_ext = ".tmp"
        else:
            file_ext = f".{file_type}"

        # Create temporary file
        temp_file = tempfile.NamedTemporaryFile(
            delete=False, suffix=file_ext, prefix=f"{panel_name}_"
        )
        temp_path = Path(temp_file.name)
        temp_file.close()

        # Download file
        logger.info(f"Downloading {panel_name} from {url}")
        response = requests.get(url, timeout=60)
        response.raise_for_status()

        with open(temp_path, "wb") as f:
            f.write(response.content)

        logger.debug(f"Downloaded {panel_name} to {temp_path}")
        return temp_path

    except Exception as e:
        logger.error(f"Error downloading {panel_name} from {url}: {e}")
        return None


def process_idtxgen_panel(file_path: Path) -> list[str]:
    """Process IDT xGen Pan-Cancer Hybridization Panel (Excel format)."""
    try:
        df = pd.read_excel(file_path)
        # Look for gene symbol column
        gene_cols = [
            col
            for col in df.columns
            if "gene" in col.lower() and "symbol" in col.lower()
        ]
        if not gene_cols:
            gene_cols = [col for col in df.columns if "gene" in col.lower()]
        if not gene_cols:
            logger.error("Could not find gene column in IDT xGen panel")
            return []

        genes = df[gene_cols[0]].dropna().astype(str).str.strip().tolist()
        genes = [g for g in genes if g and len(g) <= 20]  # Filter valid gene symbols

        return list(set(genes))  # Remove duplicates
    except Exception as e:
        logger.error(f"Error processing IDT xGen panel: {e}")
        return []


def process_illumina_trusight_panel(file_path: Path) -> list[str]:
    """Process Illumina TruSight Oncology panel (Excel format)."""
    try:
        df = pd.read_excel(file_path, skiprows=2)  # Skip first 2 rows as in R code

        # Look for gene symbol column
        gene_cols = [
            col
            for col in df.columns
            if "gene" in col.lower() and "symbol" in col.lower()
        ]
        if not gene_cols:
            logger.error("Could not find gene symbol column in TruSight panel")
            return []

        genes = df[gene_cols[0]].dropna().astype(str).str.strip()

        # Filter out specific text as in R code
        genes = genes[
            ~genes.str.contains("Small variants found in gVCF file only", na=False)
        ]

        genes_list = genes.tolist()
        genes_list = [g for g in genes_list if g and len(g) <= 20]

        return list(set(genes_list))
    except Exception as e:
        logger.error(f"Error processing Illumina TruSight panel: {e}")
        return []


def process_sureselect_allinone_panel(file_path: Path) -> list[str]:
    """Process Agilent SureSelect Cancer All-In-One Solid Tumor Assay (PDF format)."""
    try:
        reader = PdfReader(file_path)
        text = ""
        for page in reader.pages:
            text += page.extract_text() + "\n"

        lines = text.split("\n")

        # Find start and end markers as in R code
        start_idx = None
        end_idx = None

        for i, line in enumerate(lines):
            if "ABL1" in line and start_idx is None:
                start_idx = i
            if "PTPN11" in line:
                end_idx = i

        if start_idx is None or end_idx is None:
            logger.warning("Could not find gene section in SureSelect All-In-One PDF")
            return extract_genes_from_text(text)

        # Extract genes from the identified section
        gene_section = lines[start_idx : end_idx + 1]
        genes = []

        for line in gene_section:
            # Split by multiple spaces and extract potential gene symbols
            parts = re.split(r"\s{3,}", line.strip())
            for part in parts:
                if part and len(part) <= 20 and re.match(r"^[A-Z][A-Z0-9-_]*$", part):
                    genes.append(part)

        return list(set(genes))
    except Exception as e:
        logger.error(f"Error processing SureSelect All-In-One panel: {e}")
        return []


def process_sureselect_glasgow_panel(file_path: Path, panel_name: str) -> list[str]:
    """Process Agilent SureSelect Genomics Glasgow Panel (PDF format)."""
    try:
        reader = PdfReader(file_path)
        text = ""
        for page in reader.pages:
            text += page.extract_text() + "\n"

        lines = text.split("\n")

        # Determine if this is core or plus version
        is_core = "core" in panel_name.lower()
        is_plus = "plus" in panel_name.lower()

        genes = []

        if is_core:
            # Find core section (AKT1 to PMS2)
            start_idx = None
            end_idx = None

            for i, line in enumerate(lines):
                if "AKT1" in line and start_idx is None:
                    start_idx = i
                if "PMS2" in line:
                    end_idx = i

            if start_idx is not None and end_idx is not None:
                gene_section = lines[start_idx : end_idx + 1]
        elif is_plus:
            # Find plus section (ABL1 to RFX5)
            start_idx = None
            end_idx = None

            for i, line in enumerate(lines):
                if "ABL1" in line and start_idx is None:
                    start_idx = i
                if "RFX5" in line:
                    end_idx = i

            if start_idx is not None and end_idx is not None:
                gene_section = lines[start_idx : end_idx + 1]
        else:
            # Process entire document
            gene_section = lines

        for line in gene_section:
            # Clean line as in R code
            cleaned_line = re.sub(r"[*+•]", "", line)  # Remove special characters
            cleaned_line = re.sub(
                r"\s+", " ", cleaned_line
            ).strip()  # Normalize whitespace

            # Split and extract gene symbols
            parts = re.split(r"\s{2,}", cleaned_line)
            for part in parts:
                part = part.strip()
                # Apply specific corrections from R code
                if part == "TGFBRN":
                    part = "TGFBR1"

                if (
                    part
                    and part != "4"
                    and len(part) <= 20
                    and re.match(r"^[A-Z][A-Z0-9-_]*$", part)
                ):
                    genes.append(part)

        return list(set(genes))
    except Exception as e:
        logger.error(f"Error processing SureSelect Glasgow panel: {e}")
        return []


def process_ion_ampliseq_panel(file_path: Path) -> list[str]:
    """Process Ion AmpliSeq Cancer Panel (PDF format)."""
    try:
        reader = PdfReader(file_path)
        text = ""
        for page in reader.pages:
            text += page.extract_text() + "\n"

        lines = text.split("\n")

        # Find section from ABL1 to TP53
        start_idx = None
        end_idx = None

        for i, line in enumerate(lines):
            if "ABL1" in line and start_idx is None:
                start_idx = i
            if "TP53" in line:
                end_idx = i

        if start_idx is None or end_idx is None:
            logger.warning("Could not find gene section in Ion AmpliSeq PDF")
            return extract_genes_from_text(text)

        gene_section = lines[start_idx : end_idx + 1]
        genes = []

        for line in gene_section:
            # Split by multiple spaces
            parts = re.split(r"\s{2,}", line.strip())
            for part in parts:
                # Further split by single spaces to handle multiple genes per cell
                subparts = part.split()
                for subpart in subparts:
                    subpart = subpart.strip()
                    if (
                        subpart
                        and len(subpart) <= 20
                        and re.match(r"^[A-Z][A-Z0-9-_]*$", subpart)
                    ):
                        genes.append(subpart)

        return list(set(genes))
    except Exception as e:
        logger.error(f"Error processing Ion AmpliSeq panel: {e}")
        return []


def process_illumina_comprehensive_panel(file_path: Path) -> list[str]:
    """Process Illumina Comprehensive Panel v3 (ZIP with BED file format)."""
    try:
        genes = []

        if file_path.suffix.lower() == ".zip":
            # Extract ZIP file
            with zipfile.ZipFile(file_path, "r") as zip_ref:
                # Look for BED or manifest files
                bed_files = [
                    name
                    for name in zip_ref.namelist()
                    if name.endswith(".bed") or "manifest" in name.lower()
                ]

                for bed_file in bed_files:
                    with zip_ref.open(bed_file) as f:
                        content = f.read().decode("utf-8")
                        # Process as tab-delimited file
                        for line in content.split("\n"):
                            if line.strip():
                                parts = line.split("\t")
                                if len(parts) >= 4:  # BED format has at least 4 columns
                                    gene_info = parts[
                                        3
                                    ]  # 4th column typically contains gene info

                                    # Clean gene symbol as in R code
                                    gene_info = re.sub(
                                        r"OBRA_", "", gene_info
                                    )  # Remove OBRA_ prefix

                                    # Split by underscore and take first part
                                    gene = gene_info.split("_")[0]

                                    # Apply specific corrections from R code
                                    if gene == "SP":
                                        gene = "SP1"

                                    if (
                                        gene
                                        and len(gene) <= 20
                                        and re.match(r"^[A-Z][A-Z0-9-_]*$", gene)
                                    ):
                                        genes.append(gene)
        else:
            # Process as regular tab-delimited file
            with open(file_path) as f:
                for line in f:
                    if line.strip():
                        parts = line.split("\t")
                        if len(parts) >= 4:
                            gene_info = parts[3]

                            # Apply same processing as above
                            gene_info = re.sub(r"OBRA_", "", gene_info)
                            gene = gene_info.split("_")[0]

                            if gene == "SP":
                                gene = "SP1"

                            if (
                                gene
                                and len(gene) <= 20
                                and re.match(r"^[A-Z][A-Z0-9-_]*$", gene)
                            ):
                                genes.append(gene)

        return list(set(genes))
    except Exception as e:
        logger.error(f"Error processing Illumina Comprehensive panel: {e}")
        return []


def process_generic_panel(file_path: Path, file_type: str) -> list[str]:
    """Generic panel processor for unknown formats."""
    try:
        if file_type == "pdf" or file_path.suffix.lower() == ".pdf":
            reader = PdfReader(file_path)
            text = ""
            for page in reader.pages:
                text += page.extract_text() + "\n"
            return extract_genes_from_text(text)

        elif file_type in ["xlsx", "xls"] or file_path.suffix.lower() in [
            ".xlsx",
            ".xls",
        ]:
            df = pd.read_excel(file_path)
            # Look for gene-related columns
            gene_cols = [col for col in df.columns if "gene" in col.lower()]
            if gene_cols:
                genes = df[gene_cols[0]].dropna().astype(str).str.strip().tolist()
                return [g for g in genes if g and len(g) <= 20]

        return []
    except Exception as e:
        logger.error(f"Error in generic panel processing: {e}")
        return []


def extract_genes_from_text(text: str) -> list[str]:
    """Extract potential gene symbols from text using regex."""
    # Look for gene-like patterns (2-20 uppercase letters/numbers with possible hyphens)
    pattern = r"\b[A-Z][A-Z0-9-_]{1,19}\b"
    potential_genes = re.findall(pattern, text)

    # Filter common false positives
    exclude_terms = {
        "PDF",
        "PAGE",
        "FILE",
        "TABLE",
        "FIGURE",
        "PANEL",
        "TEST",
        "ANALYSIS",
        "SEQUENCING",
        "DNA",
        "RNA",
        "EXOME",
        "GENOME",
        "VARIANT",
        "MUTATION",
        "CANCER",
        "TUMOR",
        "SOLID",
        "ASSAY",
        "KIT",
        "PROTOCOL",
        "METHOD",
    }

    filtered_genes = [gene for gene in potential_genes if gene not in exclude_terms]

    return list(set(filtered_genes))


def validate_somatic_commercial_config(config: dict[str, Any]) -> list[str]:
    """
    Validate somatic commercial panels configuration.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation errors
    """
    errors = []
    somatic_config = config.get("data_sources", {}).get("somatic_commercial", {})

    if not isinstance(somatic_config, dict):
        errors.append("somatic_commercial config must be a dictionary")
        return errors

    panels_config = somatic_config.get("panels", [])
    if not isinstance(panels_config, list):
        errors.append("somatic_commercial.panels must be a list")
        return errors

    for i, panel_config in enumerate(panels_config):
        if not isinstance(panel_config, dict):
            errors.append(f"somatic_commercial.panels[{i}] must be a dictionary")
            continue

        # Check required fields
        name = panel_config.get("name")
        if not name or not isinstance(name, str):
            errors.append(
                f"somatic_commercial.panels[{i}].name is required and must be a string"
            )

        url = panel_config.get("url")
        if not url or not isinstance(url, str):
            errors.append(
                f"somatic_commercial.panels[{i}].url is required and must be a string"
            )

        # Check evidence score
        evidence_score = panel_config.get("evidence_score", 0.8)
        if not isinstance(evidence_score, int | float) or evidence_score < 0:
            errors.append(
                f"somatic_commercial.panels[{i}].evidence_score must be a non-negative number"
            )

        # Check category
        category = panel_config.get("category", "somatic")
        if category not in ["germline", "somatic"]:
            errors.append(
                f"somatic_commercial.panels[{i}].category must be either 'germline' or 'somatic'"
            )

    return errors


def get_somatic_commercial_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary information about configured somatic commercial panels.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    somatic_config = config.get("data_sources", {}).get("somatic_commercial", {})
    panels_config = somatic_config.get("panels", [])

    summary = {
        "enabled": somatic_config.get("enabled", True),
        "total_panels": len(panels_config),
        "panels": [],
        "validation_errors": validate_somatic_commercial_config(config),
    }

    for panel_config in panels_config:
        panel_summary = {
            "name": panel_config.get("name", "Unknown"),
            "url": panel_config.get("url", ""),
            "evidence_score": panel_config.get("evidence_score", 0.8),
            "category": panel_config.get("category", "somatic"),
            "type": panel_config.get("type", "auto"),
        }
        summary["panels"].append(panel_summary)

    return summary
