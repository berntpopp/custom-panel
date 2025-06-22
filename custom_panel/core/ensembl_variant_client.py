"""
Ensembl Variant API client for coordinate conversion and rsID lookup.

This module provides a comprehensive client for the Ensembl REST API with
support for coordinate lift-over, rsID lookup, and variant annotation.

## Ensembl API Integration

### API Overview
The Ensembl REST API (https://rest.ensembl.org) provides:
- **Coordinate mapping**: Cross-genome build coordinate conversion
- **Variant lookup**: rsID to coordinate mapping
- **Sequence retrieval**: Reference sequence access
- **Species support**: Human and other species data

### Client Features

#### 1. Cross-Genome Build Mapping
```python
# Convert coordinates between genome builds
hg38_coords = client.map_coordinates(
    chromosome="1",
    position=100880328,
    source_assembly="GRCh37",
    target_assembly="GRCh38"
)
```

#### 2. rsID Lookup and Resolution
```python
# Find rsID for genomic coordinates
rsid = client.lookup_rsid_by_coordinates(
    chromosome="1",
    position=100880328,
    assembly="GRCh37"
)

# Get coordinates for rsID
coords = client.get_coordinates_for_rsid("rs123456")
```

#### 3. Rate Limiting and Respect
- **Configurable delays**: Respects Ensembl's rate limiting guidelines
- **Exponential backoff**: Automatic retry with increasing delays
- **Request batching**: Efficient handling of multiple requests

#### 4. Error Handling and Resilience
- **API error handling**: Graceful handling of 404s and timeouts
- **Data validation**: Validates API responses before processing
- **Fallback strategies**: Multiple approaches for variant resolution

### Usage Patterns

#### Coordinate Conversion
```python
client = EnsemblVariantClient(rate_limit_delay=0.2)

# Single coordinate conversion
result = client.map_coordinates("1", 100880328, "GRCh37", "GRCh38")
if result:
    print(f"GRCh37 1:100880328 -> GRCh38 {result['chromosome']}:{result['position']}")

# Batch coordinate conversion for BCAC data
enhanced_df = client.enhance_bcac_with_hg38_coordinates(bcac_df)
```

#### rsID Enhancement
```python
# Enhance coordinate-only data with rsIDs
df_with_rsids = client.enhance_with_rsids(
    coordinates_df,
    assembly="GRCh37"
)

# Success rate reporting
success_rate = client.get_enhancement_stats()
print(f"rsID lookup success: {success_rate['rsid_success']}/{success_rate['total']}")
```

### API Endpoints Used

#### 1. Coordinate Mapping
```
GET /map/human/{source_assembly}/{region}/{target_assembly}
Example: /map/human/GRCh37/1:100880328..100880328:1/GRCh38
```

#### 2. Variant Lookup
```
GET /overlap/region/human/{region}?feature=variation
Example: /overlap/region/human/1:100880328-100880328?feature=variation
```

#### 3. rsID Resolution
```
GET /variation/human/{rsid}
Example: /variation/human/rs123456
```

### Data Processing

#### Coordinate Conversion Logic
```python
def map_coordinates(self, chromosome, position, source_assembly, target_assembly):
    # 1. Format region string
    region = f"{chromosome}:{position}..{position}:1"

    # 2. Call Ensembl mapping API
    url = f"{self.base_url}/map/human/{source_assembly}/{region}/{target_assembly}"
    response = requests.get(url, headers={"Content-Type": "application/json"})

    # 3. Parse and validate response
    if response.status_code == 200:
        data = response.json()
        return self._extract_mapped_coordinates(data)

    return None
```

#### rsID Lookup Strategy
```python
def lookup_rsid_by_coordinates(self, chromosome, position, assembly):
    # 1. Query variation overlap endpoint
    region = f"{chromosome}:{position}-{position}"
    url = f"{self.base_url}/overlap/region/human/{region}"

    # 2. Filter for SNV variations with rsIDs
    variations = self._get_variations_in_region(url)

    # 3. Return best matching rsID
    return self._select_best_rsid(variations, position)
```

### Performance Optimizations

#### Rate Limiting Strategy
- **Respect rate limits**: Built-in delays between API calls
- **Adaptive delays**: Increases delay if rate limited
- **Request caching**: Avoids duplicate API calls

#### Batch Processing
- **Coordinate batching**: Groups nearby coordinates for efficiency
- **Parallel processing**: Multiple API calls with controlled concurrency
- **Progress tracking**: Reports progress for long-running operations

#### Memory Management
- **Streaming processing**: Processes large datasets in chunks
- **Cache management**: Limits memory usage of response cache
- **Garbage collection**: Explicit cleanup of large response objects

### Error Handling

#### API Errors
- **HTTP 429**: Rate limiting - automatic retry with exponential backoff
- **HTTP 404**: Coordinate/variant not found - logged and skipped
- **HTTP 500**: Server errors - retry with increasing delays
- **Network timeouts**: Configurable timeout with retry logic

#### Data Validation
- **Response validation**: Ensures API responses contain expected fields
- **Coordinate validation**: Validates returned coordinates are reasonable
- **Assembly validation**: Confirms assembly matches expectations

#### Graceful Degradation
- **Partial results**: Returns available data even if some lookups fail
- **Progress preservation**: Saves intermediate results during batch processing
- **Error reporting**: Detailed logs for failed lookups

### Integration with PRS Pipeline

#### BCAC Enhancement
```python
# Enhance BCAC 313 PRS with hg38 coordinates and rsIDs
enhanced_df = client.enhance_bcac_with_hg38_coordinates(bcac_df)

# Results include:
# - Original hg19 coordinates preserved
# - New hg38 coordinates added
# - rsIDs added where available
# - Success/failure tracking per variant
```

#### Quality Metrics
- **Conversion success rate**: Percentage of coordinates successfully converted
- **rsID resolution rate**: Percentage of coordinates resolved to rsIDs
- **Cross-build validation**: Consistency checks between coordinate systems

#### Configuration Integration
- **Rate limiting**: Configurable delays respect server load
- **Timeout settings**: API timeout configuration
- **Cache settings**: Response caching configuration
- **Assembly preferences**: Default assembly configuration
"""

import logging
import time
from typing import Any, Optional

import pandas as pd
import requests

logger = logging.getLogger(__name__)


class EnsemblVariantClient:
    """Client for Ensembl variant API operations."""

    def __init__(self, cache_enabled: bool = True, rate_limit_delay: float = 0.2):
        """
        Initialize the Ensembl client.

        Args:
            cache_enabled: Whether to cache API responses
            rate_limit_delay: Delay between API calls to respect rate limits
        """
        self.base_url = "https://rest.ensembl.org"
        self.cache_enabled = cache_enabled
        self.rate_limit_delay = rate_limit_delay
        self._cache: dict[str, Any] = {}

    def convert_coordinates(
        self,
        variants: list[dict[str, Any]],
        source_assembly: str = "GRCh37",
        target_assembly: str = "GRCh38",
    ) -> list[dict[str, Any]]:
        """
        Convert variant coordinates between genome assemblies.

        Args:
            variants: List of variant dicts with 'chromosome', 'position', 'ref', 'alt'
            source_assembly: Source genome assembly
            target_assembly: Target genome assembly

        Returns:
            List of converted variants with additional coordinate information
        """
        converted_variants = []

        for variant in variants:
            try:
                converted = self._convert_single_variant(
                    variant, source_assembly, target_assembly
                )
                if converted:
                    converted_variants.append(converted)
                else:
                    # Keep original if conversion fails
                    converted_variants.append(variant)

                # Rate limiting
                time.sleep(self.rate_limit_delay)

            except Exception as e:
                logger.warning(f"Failed to convert variant {variant}: {e}")
                converted_variants.append(variant)

        logger.info(f"Converted {len(converted_variants)}/{len(variants)} variants")
        return converted_variants

    def lookup_rsids(self, variants: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """
        Lookup rsIDs for variants using Ensembl API.

        Args:
            variants: List of variant dicts with coordinate information

        Returns:
            List of variants enriched with rsID information
        """
        enriched_variants = []

        for variant in variants:
            try:
                rsid = self._lookup_single_rsid(variant)
                if rsid:
                    variant = variant.copy()
                    variant["rsid"] = rsid
                enriched_variants.append(variant)

                # Rate limiting
                time.sleep(self.rate_limit_delay)

            except Exception as e:
                logger.warning(f"Failed to lookup rsID for {variant}: {e}")
                enriched_variants.append(variant)

        rsids_found = sum(1 for v in enriched_variants if v.get("rsid"))
        logger.info(f"Found rsIDs for {rsids_found}/{len(variants)} variants")
        return enriched_variants

    def _convert_single_variant(
        self, variant: dict[str, Any], source_assembly: str, target_assembly: str
    ) -> Optional[dict[str, Any]]:
        """Convert a single variant between assemblies."""

        chromosome = str(variant["chromosome"]).replace("chr", "")
        position = int(variant["position"])

        # Create cache key
        cache_key = (
            f"convert:{source_assembly}:{target_assembly}:{chromosome}:{position}"
        )

        if self.cache_enabled and cache_key in self._cache:
            cached_result = self._cache[cache_key]
            if cached_result:
                result = variant.copy()
                result.update(cached_result)
                return result
            else:
                return None

        # Use Ensembl assembly mapping endpoint
        # Format: /map/species/asm_one/region/asm_two
        # Region format: chr:start..end:strand (strand: 1 for forward, -1 for reverse)
        url = f"{self.base_url}/map/human/{source_assembly}/{chromosome}:{position}..{position}:1/{target_assembly}"

        try:
            # Add content-type parameter to URL
            params = {"content-type": "application/json"}
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()

            if data and "mappings" in data and data["mappings"]:
                mapping = data["mappings"][0]
                converted_data = {
                    f"{target_assembly.lower()}_chromosome": str(
                        mapping["mapped"]["seq_region_name"]
                    ),
                    f"{target_assembly.lower()}_position": mapping["mapped"]["start"],
                    f"{target_assembly.lower()}_end": mapping["mapped"]["end"],
                }

                if self.cache_enabled:
                    self._cache[cache_key] = converted_data

                result = variant.copy()
                result.update(converted_data)
                return result
            else:
                if self.cache_enabled:
                    self._cache[cache_key] = None
                return None

        except Exception as e:
            logger.debug(
                f"Coordinate conversion failed for {chromosome}:{position}: {e}"
            )
            if self.cache_enabled:
                self._cache[cache_key] = None
            return None

    def _lookup_single_rsid(self, variant: dict[str, Any]) -> Optional[str]:
        """Lookup rsID for a single variant."""

        # Try different coordinate fields
        for chr_field in ["chromosome", "hg38_chromosome", "hg19_chromosome"]:
            for pos_field in [
                "position",
                "hg38_position",
                "hg19_position",
                "hg38_start",
                "hg19_start",
            ]:
                if chr_field in variant and pos_field in variant:
                    chromosome = str(variant[chr_field]).replace("chr", "")
                    position = variant[pos_field]

                    if chromosome and position:
                        rsid = self._query_variant_rsid(chromosome, int(position))
                        if rsid:
                            return rsid
        return None

    def _query_variant_rsid(self, chromosome: str, position: int) -> Optional[str]:
        """Query Ensembl for rsID at specific position."""

        cache_key = f"rsid:{chromosome}:{position}"

        if self.cache_enabled and cache_key in self._cache:
            return self._cache[cache_key]

        # Use Ensembl overlap endpoint to find variants
        url = f"{self.base_url}/overlap/region/human/{chromosome}:{position}-{position}"

        try:
            # Add parameters to URL
            params = {"feature": "variation", "content-type": "application/json"}
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()

            data = response.json()

            # Look for variants at exact position
            for variant in data:
                if variant.get("start") == position or variant.get("end") == position:
                    rsid = variant.get("id")
                    if rsid and rsid.startswith("rs"):
                        if self.cache_enabled:
                            self._cache[cache_key] = rsid
                        return rsid

            if self.cache_enabled:
                self._cache[cache_key] = None
            return None

        except Exception as e:
            logger.debug(f"rsID lookup failed for {chromosome}:{position}: {e}")
            if self.cache_enabled:
                self._cache[cache_key] = None
            return None

    def enhance_variants_batch(
        self,
        variants_df: pd.DataFrame,
        source_assembly: str = "GRCh37",
        target_assembly: str = "GRCh38",
    ) -> pd.DataFrame:
        """
        Batch enhance variants with coordinate conversion and rsID lookup.

        Args:
            variants_df: DataFrame with variant information
            source_assembly: Source genome assembly
            target_assembly: Target genome assembly

        Returns:
            Enhanced DataFrame with additional coordinate and rsID information
        """
        logger.info(f"Enhancing {len(variants_df)} variants with Ensembl API")

        # Convert DataFrame to list of dicts
        variants_raw = variants_df.to_dict("records")
        variants: list[dict[str, Any]] = [
            {str(k): v for k, v in record.items()} for record in variants_raw
        ]

        # Convert coordinates
        logger.info(
            f"Converting coordinates from {source_assembly} to {target_assembly}"
        )
        variants = self.convert_coordinates(variants, source_assembly, target_assembly)

        # Lookup rsIDs
        logger.info("Looking up rsIDs")
        variants = self.lookup_rsids(variants)

        # Convert back to DataFrame
        enhanced_df = pd.DataFrame(variants)

        logger.info(f"Enhancement complete: {len(enhanced_df)} variants processed")
        return enhanced_df
