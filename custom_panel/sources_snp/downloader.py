"""Downloader module for fetching SNP panel files from various sources."""

import hashlib
import logging
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import requests
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
from webdriver_manager.chrome import ChromeDriverManager

logger = logging.getLogger(__name__)


class PanelDownloader:
    """Handles downloading panel files from various sources."""

    def __init__(self, cache_dir: Path | None = None) -> None:
        """Initialize the downloader with optional cache directory.

        Args:
            cache_dir: Directory to cache downloaded files. If None, uses data/snp/downloads/
        """
        self.cache_dir = cache_dir or Path("data/snp/downloads")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._driver: webdriver.Chrome | None = None

    def _get_cache_filename(self, url: str, extension: str) -> Path:
        """Generate a cache filename based on URL hash.

        Args:
            url: The URL to hash
            extension: File extension to use

        Returns:
            Path to the cache file
        """
        url_hash = hashlib.md5(url.encode()).hexdigest()[:8]
        parsed = urlparse(url)
        base_name = Path(parsed.path).stem or "download"
        return self.cache_dir / f"{base_name}_{url_hash}.{extension}"

    def _download_with_requests(self, url: str, extension: str) -> tuple[Path, bool]:
        """Download file using requests library.

        Args:
            url: URL to download
            extension: File extension

        Returns:
            Tuple of (file path, was_cached)
        """
        cache_file = self._get_cache_filename(url, extension)

        if cache_file.exists():
            logger.info(f"Using cached file: {cache_file}")
            return cache_file, True

        logger.info(f"Downloading {url} to {cache_file}")

        try:
            response = requests.get(url, timeout=60, allow_redirects=True)
            response.raise_for_status()

            with open(cache_file, "wb") as f:
                f.write(response.content)

            logger.info(f"Successfully downloaded to {cache_file}")
            return cache_file, False

        except Exception as e:
            logger.error(f"Failed to download {url}: {e}")
            raise

    def _get_webdriver(self) -> webdriver.Chrome:
        """Get or create a Chrome webdriver instance."""
        if self._driver is None:
            chrome_options = Options()
            chrome_options.add_argument("--headless")
            chrome_options.add_argument("--no-sandbox")
            chrome_options.add_argument("--disable-dev-shm-usage")
            chrome_options.add_argument("--disable-gpu")
            chrome_options.add_argument("--window-size=1920,1080")
            chrome_options.add_argument("--disable-blink-features=AutomationControlled")
            chrome_options.add_experimental_option(
                "excludeSwitches", ["enable-automation"]
            )
            chrome_options.add_experimental_option("useAutomationExtension", False)

            try:
                service = Service(ChromeDriverManager().install())
                driver = webdriver.Chrome(service=service, options=chrome_options)
                self._driver = driver
                # Execute script to hide webdriver property
                driver.execute_script(
                    "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"
                )
            except Exception as e:
                logger.error(f"Failed to create webdriver: {e}")
                raise

        assert self._driver is not None
        return self._driver

    def _download_with_selenium(self, url: str, extension: str) -> tuple[Path, bool]:
        """Download file using Selenium for JavaScript-rendered pages.

        Args:
            url: URL to download
            extension: File extension

        Returns:
            Tuple of (file path, was_cached)
        """
        cache_file = self._get_cache_filename(url, extension)

        if cache_file.exists():
            logger.info(f"Using cached file: {cache_file}")
            return cache_file, True

        logger.info(f"Downloading {url} with Selenium to {cache_file}")

        try:
            driver = self._get_webdriver()
            driver.get(url)

            # Wait for page to load
            WebDriverWait(driver, 20).until(
                EC.presence_of_element_located((By.TAG_NAME, "body"))
            )

            # Additional wait for dynamic content
            time.sleep(3)

            # Get page source
            page_source = driver.page_source

            with open(cache_file, "w", encoding="utf-8") as f:
                f.write(page_source)

            logger.info(f"Successfully downloaded to {cache_file}")
            return cache_file, False

        except Exception as e:
            logger.error(f"Failed to download {url} with Selenium: {e}")
            raise

    def download(self, url: str, file_type: str, method: str = "auto") -> Path:
        """Download a file from URL using the appropriate method.

        Args:
            url: URL to download
            file_type: Type of file (pdf, html, list, tsv, etc.)
            method: Download method - "requests", "selenium", or "auto"

        Returns:
            Path to the downloaded file
        """
        if method == "auto":
            # Determine method based on file type
            if file_type == "html":
                method = "selenium"
            else:
                method = "requests"

        if method == "selenium":
            file_path, _ = self._download_with_selenium(url, file_type)
        else:
            file_path, _ = self._download_with_requests(url, file_type)

        return file_path

    def close(self) -> None:
        """Close the webdriver if it exists."""
        if self._driver:
            self._driver.quit()
            self._driver = None

    def __enter__(self) -> "PanelDownloader":
        """Context manager entry."""
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Context manager exit - close webdriver."""
        self.close()
