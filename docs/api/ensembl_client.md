# Ensembl Client

The Ensembl client provides access to genomic data from the Ensembl REST API.

## Overview

The `EnsemblClient` class handles:

- **Gene coordinate lookup** - Retrieves current genomic coordinates
- **Transcript information** - Finds canonical and MANE transcripts
- **Batch processing** - Efficiently processes multiple genes
- **Error handling** - Robust fallback mechanisms for API failures
- **Caching** - Reduces API load with intelligent caching

## Features

- **Multiple assemblies** - Support for GRCh38 and GRCh37
- **MANE transcripts** - Identifies MANE Select and Plus Clinical transcripts
- **Exon data** - Retrieves detailed exon coordinates for BED file generation
- **Rate limiting** - Respects Ensembl API guidelines
- **Retry logic** - Handles transient network issues

## API Reference

::: custom_panel.core.ensembl_client
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2