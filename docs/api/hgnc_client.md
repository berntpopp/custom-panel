# HGNC Client

The HGNC client provides gene symbol standardization and validation using the HUGO Gene Nomenclature Committee database.

## Overview

The `HGNCClient` class provides:

- **Symbol standardization** - Converts aliases to approved symbols
- **Gene validation** - Verifies gene symbols exist in HGNC
- **Batch processing** - Efficiently processes multiple genes
- **Alternative symbols** - Handles previous and alias symbols
- **Caching** - Reduces API load with local caching

## Features

- **Comprehensive search** - Searches approved symbols, aliases, and previous symbols
- **Status filtering** - Only returns approved and active genes
- **Fuzzy matching** - Handles common symbol variations
- **Error recovery** - Graceful handling of API timeouts and errors
- **Performance optimization** - Batched requests and intelligent caching

## API Reference

::: custom_panel.core.hgnc_client
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2