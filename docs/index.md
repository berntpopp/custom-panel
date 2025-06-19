# Welcome to Custom Panel!

A modern Python tool for gene panel curation and aggregation from multiple genomic databases.

## Overview

Custom Panel is a comprehensive bioinformatics tool for creating, managing, and curating gene panels for clinical genomics applications. It aggregates gene information from multiple trusted sources, applies a configurable scoring system, and generates standardized, analysis-ready output formats.

## Key Features

- **Multi-source data aggregation**: Integrates data from PanelApp, ACMG recommendations, in-house panels, manual curation lists, ClinGen, TheGenCC, COSMIC, and HPO/OMIM
- **Intelligent scoring system**: Configurable evidence weighting with veto capabilities for critical sources  
- **Gene standardization**: Automatic gene symbol standardization using HGNC
- **Genomic annotation**: Rich annotation with Ensembl coordinates, transcripts, and MANE information
- **Flexible output formats**: Excel, CSV, Parquet, and BED file generation
- **Modern architecture**: Built with Python 3.10+, Poetry, and comprehensive type hints

## Getting Started

New to Custom Panel? Start here:

- **[Installation Guide](./user_guide/installation.md)** - Set up Custom Panel on your system
- **[Running the Pipeline](./user_guide/running_pipeline.md)** - Learn the basic commands and workflows
- **[Configuration](./user_guide/configuration.md)** - Understand the scoring system and customize data sources

## Data Sources

Learn how to configure and use different data sources:

- **[COSMIC Setup](./user_guide/cosmic_setup.md)** - Configure COSMIC Cancer Gene Census access
- **[OMIM/HPO Setup](./user_guide/omim_setup.md)** - Set up OMIM and HPO data sources

## API Reference

For developers and advanced users:

- **[CLI Reference](./api/cli.md)** - Complete command-line interface documentation
- **[Core Components](./api/annotator.md)** - Gene annotation and processing engines
- **[Data Clients](./api/ensembl_client.md)** - API clients for external services

## Scientific Background

Custom Panel implements a weighted evidence aggregation system that:

1. **Source Integration**: Standardizes data from multiple sources to a common schema
2. **Gene Standardization**: Validates and standardizes gene symbols using HGNC
3. **Evidence Scoring**: Applies source-specific weights based on clinical reliability
4. **Decision Logic**: Uses configurable thresholds with veto capabilities for critical sources
5. **Genomic Annotation**: Enriches data with current coordinates and transcript information

## Support

- **Issues**: [GitHub Issues](https://github.com/berntpopp/custom-panel/issues)
- **Discussions**: [GitHub Discussions](https://github.com/berntpopp/custom-panel/discussions)
- **Source Code**: [GitHub Repository](https://github.com/berntpopp/custom-panel)