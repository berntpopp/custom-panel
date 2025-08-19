"""
COSMIC Cancer Gene Census germline data source extractor with authentication.

This module fetches and processes the COSMIC Cancer Gene Census, providing
germline evidence scoring based on tier classifications.
Includes authentication support for accessing COSMIC data through login.
"""

import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from ..core.config_manager import ConfigManager
from ..core.io import create_standard_dataframe

logger = logging.getLogger(__name__)

# COSMIC URLs
COSMIC_LOGIN_URL = "https://cancer.sanger.ac.uk/cosmic/login"
COSMIC_FILE_DOWNLOAD_URL = "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v97/cancer_gene_census.csv"


class COSMICAuthenticationError(Exception):
    """Raised when COSMIC authentication fails."""

    pass


class COSMICSession:
    """Manages COSMIC authentication and session handling."""

    def __init__(self, email: str, password: str):
        """
        Initialize COSMIC session.

        Args:
            email: COSMIC account email
            password: COSMIC account password
        """
        self.email = email
        self.password = password
        self.session = requests.Session()
        self._setup_session()

    def _setup_session(self) -> None:
        """Setup session with appropriate headers."""
        self.session.headers.update(
            {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/137.0.0.0 Safari/537.36",
                "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
                "Accept-Language": "en-US,en;q=0.9",
                "Accept-Encoding": "gzip, deflate, br",
                "Cache-Control": "max-age=0",
                "Sec-CH-UA": '"Google Chrome";v="137", "Chromium";v="137", "Not/A)Brand";v="24"',
                "Sec-CH-UA-Mobile": "?0",
                "Sec-CH-UA-Platform": '"Windows"',
                "Sec-Fetch-Dest": "document",
                "Sec-Fetch-Mode": "navigate",
                "Sec-Fetch-Site": "same-origin",
                "Sec-Fetch-User": "?1",
                "Upgrade-Insecure-Requests": "1",
            }
        )

    def login(self) -> None:
        """
        Authenticate with COSMIC using credentials.

        Raises:
            COSMICAuthenticationError: If login fails
        """
        logger.info("Authenticating with COSMIC...")

        try:
            # First, get the login page to establish session
            login_page_response = self.session.get(COSMIC_LOGIN_URL, timeout=30)
            login_page_response.raise_for_status()

            # Prepare login data
            login_data = {
                "email": self.email,
                "pass": self.password,
                "r_url": "",  # redirect URL (empty for default)
                "d": "0",  # download flag
            }

            # Set appropriate headers for form submission
            login_headers = {
                "Content-Type": "application/x-www-form-urlencoded",
                "Origin": "https://cancer.sanger.ac.uk",
                "Referer": COSMIC_LOGIN_URL,
            }

            # Perform login
            login_response = self.session.post(
                COSMIC_LOGIN_URL,
                data=login_data,
                headers=login_headers,
                timeout=30,
                allow_redirects=True,
            )
            login_response.raise_for_status()

            # Check if login was successful
            logger.debug(f"Login response URL: {login_response.url}")
            logger.debug(f"Login response status: {login_response.status_code}")

            # Check for explicit error indicators
            response_text_lower = login_response.text.lower()
            if "invalid" in response_text_lower or "incorrect" in response_text_lower:
                raise COSMICAuthenticationError("Invalid credentials")

            # If we're still on the login page, check for error messages
            if "login" in login_response.url.lower():
                if "error" in response_text_lower or "failed" in response_text_lower:
                    raise COSMICAuthenticationError(
                        "Login failed - please check credentials"
                    )
                # Sometimes COSMIC returns to login page even on success
                logger.warning("Still on login page, but no explicit error detected")

            logger.info("Successfully authenticated with COSMIC")

        except requests.RequestException as e:
            logger.error(f"Network error during COSMIC login: {e}")
            raise COSMICAuthenticationError(f"Failed to connect to COSMIC: {e}") from e
        except Exception as e:
            logger.error(f"Unexpected error during COSMIC login: {e}")
            raise COSMICAuthenticationError(f"Login failed: {e}") from e

    def get_download_url(self, file_endpoint: str) -> str:
        """
        Get authenticated download URL for a COSMIC file.

        Args:
            file_endpoint: The file download endpoint

        Returns:
            Authenticated download URL

        Raises:
            COSMICAuthenticationError: If unable to get download URL
        """
        try:
            logger.info(f"Getting download URL for: {file_endpoint}")

            response = self.session.get(file_endpoint, timeout=30)
            response.raise_for_status()

            logger.debug(f"Download URL response status: {response.status_code}")
            logger.debug(
                f"Download URL response content type: {response.headers.get('content-type', 'unknown')}"
            )

            # Check if we're being redirected to login (authentication failed)
            if "login" in response.url.lower():
                raise COSMICAuthenticationError(
                    "Authentication expired or invalid - redirected to login"
                )

            # Parse JSON response to get the actual download URL
            try:
                data = response.json()
                download_url = data.get("url")
                if not download_url:
                    logger.error(f"No download URL in response: {data}")
                    raise COSMICAuthenticationError(
                        "No download URL returned by COSMIC"
                    )

                logger.info("Successfully obtained authenticated download URL")
                logger.debug(
                    f"Download URL: {download_url[:100]}..."
                )  # Log first 100 chars for security
                return download_url

            except ValueError as e:
                logger.error(
                    f"Invalid JSON response from COSMIC. Response text (first 500 chars): {response.text[:500]}"
                )
                logger.error(f"Full response headers: {dict(response.headers)}")
                raise COSMICAuthenticationError(f"Invalid response format: {e}") from e

        except requests.RequestException as e:
            logger.error(f"Network error getting download URL: {e}")
            raise COSMICAuthenticationError(f"Failed to get download URL: {e}") from e

    def download_file(self, download_url: str, cache_path: Path) -> None:
        """
        Download file from authenticated URL.

        Args:
            download_url: Authenticated download URL
            cache_path: Local path to save the file

        Raises:
            COSMICAuthenticationError: If download fails
        """
        logger.info(f"Downloading COSMIC file to: {cache_path}")

        # Create cache directory
        cache_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            response = self.session.get(download_url, timeout=300, stream=True)
            response.raise_for_status()

            # Save file with progress tracking
            total_size = int(response.headers.get("content-length", 0))
            downloaded_size = 0

            with open(cache_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded_size += len(chunk)

                        # Log progress for large files
                        if (
                            total_size > 0 and downloaded_size % (1024 * 1024) == 0
                        ):  # Every MB
                            progress = (downloaded_size / total_size) * 100
                            logger.debug(f"Download progress: {progress:.1f}%")

            logger.info(
                f"Successfully downloaded COSMIC file ({downloaded_size:,} bytes)"
            )

        except requests.RequestException as e:
            logger.error(f"Failed to download COSMIC file: {e}")
            if cache_path.exists():
                cache_path.unlink()  # Remove partial download
            raise COSMICAuthenticationError(f"Download failed: {e}") from e


def _download_cosmic_census_authenticated(
    email: str, password: str, cache_path: Path
) -> None:
    """
    Download COSMIC Cancer Gene Census file with authentication.

    Args:
        email: COSMIC account email
        password: COSMIC account password
        cache_path: Local path to save the file

    Raises:
        COSMICAuthenticationError: If authentication or download fails
    """
    try:
        # Create authenticated session
        cosmic_session = COSMICSession(email, password)

        # Login to COSMIC
        cosmic_session.login()

        # Get authenticated download URL
        download_url = cosmic_session.get_download_url(COSMIC_FILE_DOWNLOAD_URL)

        # Download the file
        cosmic_session.download_file(download_url, cache_path)

    except COSMICAuthenticationError:
        raise
    except Exception as e:
        logger.error(f"Unexpected error during COSMIC download: {e}")
        raise COSMICAuthenticationError(f"Download failed: {e}") from e


def _download_cosmic_census(url: str, cache_path: Path) -> None:
    """
    Download COSMIC Cancer Gene Census file (legacy unauthenticated method).

    Args:
        url: URL to download from
        cache_path: Local path to save the file

    Raises:
        requests.RequestException: If download fails
    """
    logger.info(f"Downloading COSMIC Cancer Gene Census from: {url}")

    # Create cache directory
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    # Download with proper headers
    headers = {
        "User-Agent": "custom-panel/1.0 (Python scientific tool for gene panel curation)"
    }

    try:
        response = requests.get(url, headers=headers, timeout=60, stream=True)
        response.raise_for_status()

        # Save file
        with open(cache_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        logger.info(f"Successfully downloaded COSMIC census to: {cache_path}")

    except requests.RequestException as e:
        logger.error(f"Failed to download COSMIC census: {e}")
        raise


def _load_cosmic_census(cache_path: Path) -> pd.DataFrame:
    """
    Load and validate COSMIC Cancer Gene Census file (CSV or TSV format).

    Args:
        cache_path: Path to the cached CSV or TSV file

    Returns:
        Loaded DataFrame

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid
    """
    if not cache_path.exists():
        raise FileNotFoundError(f"COSMIC census file not found: {cache_path}")

    try:
        # Determine separator based on file extension
        if cache_path.suffix.lower() == ".tsv":
            separator = "\t"
            logger.info(f"Loading COSMIC census TSV file: {cache_path}")
        else:
            separator = ","
            logger.info(f"Loading COSMIC census CSV file: {cache_path}")

        df = pd.read_csv(cache_path, sep=separator)
        logger.info(f"Loaded COSMIC census with {len(df)} genes")

        # Validate required columns
        required_columns = ["Gene Symbol", "Tier"]
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(
                f"COSMIC census missing required columns: {missing_columns}"
            )

        # Check for germline and somatic columns (may vary by COSMIC version)
        # Try different possible column names
        germline_cols = ["Germline", "Germline Mutation"]
        somatic_cols = ["Somatic", "Somatic Mutation", "Tumour Types(Somatic)"]

        germline_col = None
        somatic_col = None

        for col in germline_cols:
            if col in df.columns:
                germline_col = col
                break

        for col in somatic_cols:
            if col in df.columns:
                somatic_col = col
                break

        if not germline_col and not somatic_col:
            logger.warning("No germline or somatic columns found in COSMIC census")
        else:
            logger.info(
                f"Found germline column: {germline_col}, somatic column: {somatic_col}"
            )

        # Store column names for later use
        df.attrs["germline_col"] = germline_col
        df.attrs["somatic_col"] = somatic_col

        return df

    except Exception as e:
        logger.error(f"Failed to load COSMIC census: {e}")
        raise ValueError(f"Invalid COSMIC census format: {e}") from e


def _find_cosmic_census_file(cache_dir: Path) -> Path | None:
    """
    Find COSMIC census file in cache directory (CSV or TSV format).

    Args:
        cache_dir: Cache directory to search

    Returns:
        Path to found file or None if not found
    """
    # Check for both CSV and TSV files
    csv_path = cache_dir / "cosmic_gene_census.csv"
    tsv_path = cache_dir / "cosmic_gene_census.tsv"

    if tsv_path.exists():
        return tsv_path
    elif csv_path.exists():
        return csv_path
    else:
        return None


def _is_cache_valid(cache_path: Path, expiry_days: int) -> bool:
    """
    Check if cached file is still valid.

    Args:
        cache_path: Path to cached file
        expiry_days: Number of days before cache expires

    Returns:
        True if cache is valid, False otherwise
    """
    if not cache_path.exists():
        return False

    file_age = datetime.now() - datetime.fromtimestamp(cache_path.stat().st_mtime)
    return file_age < timedelta(days=expiry_days)


def _calculate_cosmic_score(tier: str, tier_weights: dict[str, float]) -> float:
    """
    Calculate evidence score based on COSMIC tier.

    Args:
        tier: COSMIC tier classification
        tier_weights: Mapping of tiers to weights

    Returns:
        Evidence score (0.0-1.0)
    """
    # Clean tier string
    tier = str(tier).strip() if pd.notna(tier) else ""

    # Get weight, defaulting to unknown tier weight
    return tier_weights.get(tier, tier_weights.get("", 0.4))


def _process_cosmic_genes(
    df: pd.DataFrame, category: str, category_config: dict[str, Any]
) -> pd.DataFrame:
    """
    Process COSMIC genes for a specific category (germline or somatic).

    Args:
        df: COSMIC census DataFrame
        category: Either "germline" or "somatic"
        category_config: Configuration for this category

    Returns:
        Standardized DataFrame for this category
    """
    if not category_config.get("enabled", False):
        logger.info(f"COSMIC {category} scoring is disabled")
        return pd.DataFrame()

    # Get the appropriate column name
    col_name = df.attrs.get(f"{category}_col")
    if not col_name:
        logger.warning(f"No {category} column found in COSMIC census")
        return pd.DataFrame()

    # Filter for genes relevant to this category
    # Look for "yes", "y", or non-empty values indicating presence
    category_df = (
        df[df[col_name].notna()].copy() if col_name in df.columns else df.copy()
    )
    category_df = category_df[
        category_df[col_name].astype(str).str.lower().isin(["yes", "y"])
        | (category_df[col_name].astype(str).str.len() > 0)
    ]

    if category_df.empty:
        logger.warning(f"No {category} genes found in COSMIC census")
        return pd.DataFrame()

    logger.info(f"Found {len(category_df)} COSMIC {category} genes")

    # Calculate evidence scores based on tier weights
    tier_weights = category_config.get("tier_weights", {"": 0.4})
    evidence_scores = [
        _calculate_cosmic_score(tier, tier_weights) for tier in category_df["Tier"]
    ]

    # Create source details
    source_details = [
        f"Tier:{tier}|Category:{category}|Date:{datetime.now().strftime('%Y-%m-%d')}"
        for tier in category_df["Tier"]
    ]

    # Create standardized DataFrame
    genes = category_df["Gene Symbol"].tolist()
    source_name = f"COSMIC_{category.title()}"

    result_df = create_standard_dataframe(
        genes=genes,
        source_name=source_name,
        evidence_scores=evidence_scores,
        source_details=source_details,
        gene_names_reported=genes,
    )

    return result_df


def fetch_cosmic_germline_data(config: dict[str, Any]) -> pd.DataFrame:
    """
    Fetch COSMIC Cancer Gene Census germline data with authentication and caching.
    Focuses exclusively on germline variants.

    Args:
        config: Configuration dictionary

    Returns:
        Standardized DataFrame with COSMIC germline data only
    """
    config_manager = ConfigManager(config)
    cosmic_config = config_manager.get_source_config("COSMIC_Germline")

    if not cosmic_config.get("enabled", False):
        logger.info("COSMIC data source is disabled")
        return pd.DataFrame()

    cache_dir = Path(cosmic_config.get("cache_dir", ".cache/cosmic"))
    cache_expiry_days = cosmic_config.get("cache_expiry_days", 30)

    # Look for existing cache file (CSV or TSV)
    cached_file = _find_cosmic_census_file(cache_dir)

    if cached_file and _is_cache_valid(cached_file, cache_expiry_days):
        logger.info(f"Using valid cached COSMIC file: {cached_file}")
        cache_path = cached_file
    else:
        # Default to CSV for downloads
        cache_path = cache_dir / "cosmic_gene_census.csv"
        logger.info("COSMIC cache expired or missing, downloading new data...")

        # Try authenticated download first
        email = cosmic_config.get("email")
        password = cosmic_config.get("password")

        if email and password:
            logger.info("Using authenticated COSMIC download")
            try:
                _download_cosmic_census_authenticated(email, password, cache_path)
            except COSMICAuthenticationError as e:
                logger.error(f"Authenticated download failed: {e}")
                logger.error("Please check your COSMIC credentials in config.local.yml")
                logger.error("Required configuration:")
                logger.error("data_sources:")
                logger.error("  cosmic:")
                logger.error("    email: your-cosmic-email@example.com")
                logger.error("    password: your-cosmic-password")
                return pd.DataFrame()
        else:
            # Fallback to legacy URL if provided
            census_url = cosmic_config.get("census_url")
            if census_url:
                logger.warning(
                    "No COSMIC credentials found, trying legacy URL download"
                )
                try:
                    _download_cosmic_census(census_url, cache_path)
                except requests.RequestException as e:
                    logger.error(f"Legacy download failed: {e}")
                    # Check for any existing cache file (even if expired)
                    cached_file = _find_cosmic_census_file(cache_dir)
                    if cached_file:
                        logger.warning(f"Using expired cache file: {cached_file}")
                        cache_path = cached_file
                    else:
                        logger.error("No cache file available, cannot proceed")
                        return pd.DataFrame()
            else:
                logger.error("No COSMIC credentials or legacy URL configured")
                # Check for any existing cache file as final fallback
                cached_file = _find_cosmic_census_file(cache_dir)
                if cached_file:
                    logger.warning(
                        f"Using existing cache file despite configuration issues: {cached_file}"
                    )
                    cache_path = cached_file
                else:
                    logger.error("Please add COSMIC credentials to config.local.yml:")
                    logger.error("data_sources:")
                    logger.error("  cosmic:")
                    logger.error("    enabled: true")
                    logger.error("    email: your-cosmic-email@example.com")
                    logger.error("    password: your-cosmic-password")
                    return pd.DataFrame()

    # Load the census data
    try:
        df = _load_cosmic_census(cache_path)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Failed to load COSMIC census: {e}")
        return pd.DataFrame()

    # Process germline category only
    germline_config = cosmic_config.get("germline_scoring", {})

    # Ensure germline scoring is enabled (default to true if section exists)
    if not germline_config.get("enabled", True):
        logger.warning(
            "COSMIC germline scoring is disabled. Enable it in configuration to include COSMIC data."
        )
        return pd.DataFrame()

    germline_df = _process_cosmic_genes(df, "germline", germline_config)

    # Add a category column for the merger
    if not germline_df.empty:
        germline_df["category"] = "germline"

    if germline_df.empty:
        logger.warning("No COSMIC germline data processed")
        return pd.DataFrame()

    logger.info(f"Created COSMIC germline dataset with {len(germline_df)} gene records")
    return germline_df


def validate_cosmic_config(config: dict[str, Any]) -> list[str]:
    """
    Validate COSMIC germline configuration.

    Args:
        config: Configuration dictionary

    Returns:
        List of validation errors
    """
    errors: list[str] = []
    config_manager = ConfigManager(config)
    cosmic_config = config_manager.get_source_config("COSMIC_Germline")

    if not cosmic_config.get("enabled", False):
        return errors  # Skip validation if disabled

    # Check authentication credentials
    email = cosmic_config.get("email")
    password = cosmic_config.get("password")
    census_url = cosmic_config.get("census_url")

    if not email and not password and not census_url:
        errors.append(
            "COSMIC requires either credentials (email/password) or legacy census_url"
        )
    elif email and not password:
        errors.append("COSMIC email provided but password missing")
    elif password and not email:
        errors.append("COSMIC password provided but email missing")
    elif email and "@" not in email:
        errors.append("COSMIC email appears to be invalid")

    # Validate cache settings
    cache_expiry = cosmic_config.get("cache_expiry_days")
    if cache_expiry is not None and (
        not isinstance(cache_expiry, int) or cache_expiry < 1
    ):
        errors.append("COSMIC cache_expiry_days must be a positive integer")

    # Validate scoring configurations
    category_config = cosmic_config.get("germline_scoring", {})
    if category_config.get("enabled", True):
        tier_weights = category_config.get("tier_weights", {})
        if not isinstance(tier_weights, dict):
            errors.append("COSMIC germline_scoring tier_weights must be a dictionary")
        else:
            for tier, weight in tier_weights.items():
                if not isinstance(weight, int | float) or weight < 0 or weight > 1:
                    errors.append(
                        f"COSMIC germline_scoring tier weight for '{tier}' must be between 0 and 1"
                    )

    return errors


def get_cosmic_summary(config: dict[str, Any]) -> dict[str, Any]:
    """
    Get summary of COSMIC germline configuration and cached data.

    Args:
        config: Configuration dictionary

    Returns:
        Summary dictionary
    """
    config_manager = ConfigManager(config)
    cosmic_config = config_manager.get_source_config("COSMIC_Germline")

    summary = {
        "enabled": cosmic_config.get("enabled", False),
        "has_credentials": bool(
            cosmic_config.get("email") and cosmic_config.get("password")
        ),
        "census_url": cosmic_config.get("census_url"),
        "cache_dir": cosmic_config.get("cache_dir", ".cache/cosmic"),
        "cache_expiry_days": cosmic_config.get("cache_expiry_days", 30),
        "germline_enabled": config_manager.get_nested(
            "data_sources",
            "COSMIC_Germline",
            "germline_scoring",
            "enabled",
            default=True,
        ),
        "validation_errors": validate_cosmic_config(config),
    }

    # Check cache status
    cache_path = Path(summary["cache_dir"]) / "cosmic_gene_census.csv"
    summary["cache_exists"] = cache_path.exists()
    if cache_path.exists():
        cache_age = datetime.now() - datetime.fromtimestamp(cache_path.stat().st_mtime)
        summary["cache_age_days"] = cache_age.days
        summary["cache_valid"] = cache_age < timedelta(
            days=summary["cache_expiry_days"]
        )
    else:
        summary["cache_age_days"] = None
        summary["cache_valid"] = False

    return summary
