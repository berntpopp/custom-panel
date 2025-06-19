# Gene Annotation Engine

The annotation engine is responsible for enriching gene panels with genomic coordinates, transcript information, and gene descriptions.

## Overview

The `GeneAnnotator` class provides comprehensive gene annotation functionality:

- **HGNC standardization** - Validates and standardizes gene symbols
- **Genomic coordinates** - Retrieves current coordinates from Ensembl
- **Transcript information** - Identifies canonical and MANE transcripts
- **Gene descriptions** - Adds descriptive information for each gene
- **Quality control** - Validates data consistency and completeness

## API Reference

::: custom_panel.engine.annotator
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2