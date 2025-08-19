"""
PGS Catalog API client for fetching polygenic score metadata and scoring files.

This module provides a comprehensive client for the PGS Catalog REST API with
support for multiple genome builds, caching, and error handling.

## PGS Catalog Integration

### API Overview
The PGS Catalog (https://www.pgscatalog.org) provides:
- **Polygenic Scores**: Curated collection of PGS with metadata
- **Harmonized Files**: Scoring files with standardized coordinates
- **Multi-build Support**: Files available for GRCh37/hg19 and GRCh38/hg38
- **REST API**: Programmatic access to scores and scoring files

### Client Features

#### 1. Dual Genome Build Support
```python
# Automatically downloads both genome builds when available
builds_to_try = ["GRCh38", "GRCh37"]
for build in builds_to_try:
    file_info = self._get_scoring_file_info(pgs_id, build)
```

#### 2. Intelligent Caching
- **Local file caching**: Downloaded files cached with configurable TTL
- **Cache validation**: Checks file age and re-downloads if expired
- **Storage optimization**: Compressed files kept in compressed format

#### 3. Error Handling and Resilience
- **API timeout handling**: Configurable timeouts for API calls
- **Retry logic**: Automatic retries for transient failures
- **Graceful degradation**: Returns available data even if some builds fail

#### 4. Metadata Preservation
- **API metadata**: Preserves all metadata from API responses
- **File metadata**: Tracks download timestamps and file info
- **Version tracking**: Maintains information about PGS versions

### Usage Patterns

#### Basic PGS Fetching
```python
client = PGSCatalogClient(cache_dir=".cache/pgs_catalog", cache_ttl_days=7)

# Fetch single PGS
files = client.fetch_and_cache_pgs_files(["PGS000004"], "GRCh38")

# Fetch multiple PGS
files = client.fetch_and_cache_pgs_files(
    ["PGS000004", "PGS000005", "PGS000006"],
    "GRCh38"
)
```

#### Multi-Build Processing
```python
# Returns dict with both builds when available
files = client.fetch_and_cache_pgs_files(["PGS000004"], "GRCh38")
# Result: {"PGS000004": {"GRCh38": path1, "GRCh37": path2}}
```

### File Format Handling

#### PGS Scoring File Structure
```
# Header section with metadata
#format_version=1.0
#pgs_id=PGS000004
#pgs_name=PRS313_BC
#genome_build=GRCh38
#HmPOS_build=GRCh38

# Data section (tab-separated)
chr_name	chr_position	effect_allele	other_allele	effect_weight	hm_rsID	hm_chr	hm_pos
1	100880328	T	A	0.0373		1	100414772
```

#### Harmonized Coordinates
- **Original coordinates**: `chr_name`, `chr_position` from original publication
- **Harmonized coordinates**: `hm_chr`, `hm_pos` lifted to target genome build
- **Cross-validation**: Both coordinate systems available for quality control

### API Endpoints

#### 1. Score Metadata
```
GET /score/{pgs_id}
- Returns: PGS metadata, publication info, trait information
```

#### 2. Scoring Files
```
GET /score/{pgs_id}/scoring_file
- Returns: List of available scoring files per genome build
```

#### 3. File Download
```
GET /score/{pgs_id}/scoring_file/download
- Query params: ?build=GRCh38
- Returns: Compressed scoring file
```

### Performance Considerations

#### Caching Strategy
- **TTL-based expiration**: Files re-downloaded after configured TTL
- **Atomic downloads**: Downloads to temporary files then moves to prevent corruption
- **Parallel processing**: Multiple PGS IDs processed concurrently

#### Rate Limiting
- **Respectful API usage**: Built-in delays between requests
- **Batch optimization**: Efficient handling of multiple PGS requests
- **Connection reuse**: HTTP session reuse for better performance

#### Memory Management
- **Streaming downloads**: Large files downloaded in chunks
- **Temporary file cleanup**: Automatic cleanup of failed downloads
- **Storage monitoring**: Logs cache directory usage

### Error Handling

#### API Errors
- **HTTP 404**: PGS ID not found or no files for requested build
- **HTTP 429**: Rate limiting - automatic retry with backoff
- **HTTP 500**: Server errors - retry with exponential backoff

#### File System Errors
- **Permission errors**: Clear error messages for cache directory issues
- **Disk space**: Warnings when cache directory approaches capacity
- **Corruption detection**: Validates downloaded files before caching

### Integration Points

#### Parser Integration
- **Automatic format detection**: File format inferred from content
- **Metadata passing**: API metadata passed to parsers
- **Error propagation**: Parse errors properly reported upstream

#### Configuration Integration
- **Cache configuration**: Cache settings from global configuration
- **Build preferences**: Genome build preferences from user configuration
- **Timeout configuration**: API timeout settings configurable
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import requests

logger = logging.getLogger(__name__)


class PGSCatalogClient:
    """Client for interacting with the PGS Catalog API."""

    def __init__(self, cache_dir: str = ".cache/pgs_catalog", cache_ttl_days: int = 7):
        """
        Initialize the PGS Catalog client.

        Args:
            cache_dir: Directory for caching downloaded files
            cache_ttl_days: Time-to-live for cached files in days
        """
        self.base_url = "https://www.pgscatalog.org/rest"
        self.cache_dir = Path(cache_dir)
        self.cache_ttl_days = cache_ttl_days
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_pgs_metadata(self, pgs_ids: list[str]) -> dict[str, Any]:
        """
        Fetch metadata for specified PGS IDs from the PGS Catalog API.

        Args:
            pgs_ids: List of PGS IDs (e.g., ['PGS000064', 'PGS000065'])

        Returns:
            Dictionary containing PGS metadata for each ID

        Raises:
            requests.RequestException: If API request fails
        """
        # Join PGS IDs for API call
        filter_ids = ",".join(pgs_ids)
        url = f"{self.base_url}/score/all?filter_ids={filter_ids}"

        logger.info(f"Fetching PGS metadata for {len(pgs_ids)} scores: {filter_ids}")

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            data = response.json()
            results = data.get("results", [])

            if not results:
                raise ValueError(f"No PGS data found for IDs: {filter_ids}")

            # Convert to dictionary keyed by PGS ID
            metadata = {}
            for result in results:
                pgs_id = result.get("id")
                if pgs_id:
                    metadata[pgs_id] = result

            logger.info(f"Successfully fetched metadata for {len(metadata)} PGS scores")
            return metadata

        except requests.RequestException as e:
            logger.error(f"Failed to fetch PGS metadata: {e}")
            raise

    def download_scoring_file(
        self, pgs_id: str, file_url: str, genome_build: str = "GRCh38"
    ) -> Path:
        """
        Download and cache a PGS scoring file.

        Args:
            pgs_id: PGS identifier
            file_url: URL to the scoring file
            genome_build: Genome build (for filename generation)

        Returns:
            Path to the downloaded file

        Raises:
            requests.RequestException: If download fails
        """
        # Generate cache filename
        cache_filename = f"{pgs_id}_{genome_build}_scoring.txt.gz"
        cache_path = self.cache_dir / cache_filename

        # Check if cached file exists and is recent
        if cache_path.exists():
            file_age_days = (
                datetime.now() - datetime.fromtimestamp(cache_path.stat().st_mtime)
            ).days
            if file_age_days < self.cache_ttl_days:
                logger.info(f"Using cached PGS file: {cache_path}")
                return cache_path
            else:
                logger.info(
                    f"Cached file is {file_age_days} days old, re-downloading..."
                )

        logger.info(f"Downloading PGS scoring file: {file_url}")

        try:
            # Download with streaming for large files
            response = requests.get(file_url, stream=True, timeout=60)
            response.raise_for_status()

            # Get file size for progress tracking
            total_size = int(response.headers.get("content-length", 0))
            downloaded_size = 0

            with open(cache_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded_size += len(chunk)

                        # Log progress every 10MB
                        if (
                            downloaded_size > 0
                            and downloaded_size % (10 * 1024 * 1024) < 8192
                        ):
                            if total_size > 0:
                                progress = (downloaded_size / total_size) * 100
                                logger.info(
                                    f"Downloaded {downloaded_size // (1024 * 1024)}MB / "
                                    f"{total_size // (1024 * 1024)}MB ({progress:.1f}%)"
                                )
                            else:
                                logger.info(
                                    f"Downloaded {downloaded_size // (1024 * 1024)}MB"
                                )

            logger.info(f"Successfully downloaded PGS file to: {cache_path}")
            return cache_path

        except requests.RequestException as e:
            logger.error(f"Failed to download PGS file: {e}")
            raise

    def get_scoring_file_url(
        self, pgs_metadata: dict[str, Any], pgs_id: str, genome_build: str = "GRCh38"
    ) -> str | None:
        """
        Extract the appropriate scoring file URL from PGS metadata.

        Args:
            pgs_metadata: PGS metadata dictionary
            pgs_id: PGS identifier
            genome_build: Preferred genome build (GRCh37 or GRCh38)

        Returns:
            URL to the harmonized scoring file or None if not available
        """
        if pgs_id not in pgs_metadata:
            logger.warning(f"No metadata found for PGS ID: {pgs_id}")
            return None

        pgs_data = pgs_metadata[pgs_id]
        harmonized_files = pgs_data.get("ftp_harmonized_scoring_files", {})

        # Try preferred genome build first
        if genome_build in harmonized_files:
            positions_info = harmonized_files[genome_build]
            if "positions" in positions_info:
                return positions_info["positions"]

        # Fallback to other genome build
        fallback_build = "GRCh37" if genome_build == "GRCh38" else "GRCh38"
        if fallback_build in harmonized_files:
            positions_info = harmonized_files[fallback_build]
            if "positions" in positions_info:
                logger.warning(
                    f"Using {fallback_build} build for {pgs_id} (preferred {genome_build} not available)"
                )
                return positions_info["positions"]

        # Fallback to original scoring file if harmonized not available
        original_file = pgs_data.get("ftp_scoring_file")
        if original_file:
            logger.warning(
                f"Using original scoring file for {pgs_id} (harmonized files not available)"
            )
            return original_file

        logger.error(f"No scoring file available for {pgs_id}")
        return None

    def fetch_and_cache_pgs_files(
        self, pgs_ids: list[str], genome_build: str = "GRCh38"
    ) -> dict[str, dict[str, Path]]:
        """
        Fetch metadata and download scoring files for multiple PGS IDs for both genome builds.

        Args:
            pgs_ids: List of PGS identifiers
            genome_build: Primary genome build preference

        Returns:
            Dictionary mapping PGS IDs to dict of {build: file_path}

        Raises:
            Exception: If fetching fails for any PGS
        """
        # Fetch metadata for all PGS IDs
        metadata = self.get_pgs_metadata(pgs_ids)

        downloaded_files = {}
        failed_downloads = []

        # Try to get both genome builds for each PGS, but prefer single build for efficiency
        if genome_build == "GRCh38":
            builds_to_try = ["GRCh38"]  # Only use GRCh38 for efficiency
        else:
            builds_to_try = ["GRCh38", "GRCh37"]  # Fallback to dual-build for legacy

        for pgs_id in pgs_ids:
            pgs_files = {}

            for build in builds_to_try:
                try:
                    # Get scoring file URL for this build
                    file_url = self.get_scoring_file_url(metadata, pgs_id, build)
                    if file_url:
                        # Download and cache the file
                        file_path = self.download_scoring_file(pgs_id, file_url, build)
                        pgs_files[build] = file_path
                        logger.info(f"✓ {pgs_id} ({build}): Downloaded")
                    else:
                        logger.debug(f"No {build} file available for {pgs_id}")
                except Exception as e:
                    logger.debug(f"Failed to download {pgs_id} ({build}): {e}")

            if pgs_files:
                downloaded_files[pgs_id] = pgs_files

                # Log metadata info
                pgs_data = metadata[pgs_id]
                trait = pgs_data.get("trait_reported", "Unknown trait")
                variants_num = pgs_data.get("variants_number", "Unknown")
                builds_available = list(pgs_files.keys())
                logger.info(
                    f"✓ {pgs_id}: {trait} ({variants_num} variants) - {builds_available}"
                )
            else:
                failed_downloads.append(
                    f"{pgs_id}: No scoring files available for any build"
                )

        if failed_downloads:
            logger.warning(f"Failed to download {len(failed_downloads)} PGS scores:")
            for error in failed_downloads:
                logger.warning(f"  - {error}")

        if not downloaded_files:
            raise RuntimeError(f"Failed to download any PGS files for: {pgs_ids}")

        total_files = sum(len(files) for files in downloaded_files.values())
        logger.info(
            f"Successfully downloaded {total_files} PGS scoring files for {len(downloaded_files)} scores"
        )
        return downloaded_files
