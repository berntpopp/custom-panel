# Custom Panel

[![CI](https://github.com/berntpopp/custom-panel/actions/workflows/docs.yml/badge.svg)](https://github.com/berntpopp/custom-panel/actions)
[![Documentation](https://img.shields.io/badge/docs-view%20site-blue.svg)](https://berntpopp.github.io/custom-panel/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A modern Python tool for gene panel curation and aggregation from multiple genomic databases.

`custom-panel` is a comprehensive bioinformatics tool for creating, managing, and curating gene panels. It aggregates gene information from trusted sources, applies a configurable scoring system, and generates standardized, analysis-ready output formats.

## Key Features

- **Multi-source data aggregation**: Integrates data from PanelApp, ACMG recommendations, in-house panels, manual curation lists, ClinGen, TheGenCC, COSMIC, and HPO/OMIM
- **Intelligent scoring system**: Configurable evidence weighting with veto capabilities for critical sources
- **Gene standardization**: Automatic gene symbol standardization using HGNC
- **Genomic annotation**: Rich annotation with Ensembl coordinates, transcripts, and MANE information
- **Flexible output formats**: Excel, CSV, Parquet, and BED file generation
- **Modern architecture**: Built with Python 3.10+, Poetry, and comprehensive type hints

## Quick Start

```bash
# Install with Poetry
poetry install

# Run the complete pipeline
custom-panel run --output-dir results

# Check configuration
custom-panel config-check
```

## Documentation

**📚 For detailed installation, usage, configuration, and API references, please visit the [full documentation site](https://berntpopp.github.io/custom-panel/).**

The documentation includes:
- Complete installation and setup guides
- Configuration examples and best practices
- Data source setup (COSMIC, OMIM/HPO)
- API reference for developers
- Scientific methodology and quality control details

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/berntpopp/custom-panel/issues)
- **Discussions**: [GitHub Discussions](https://github.com/berntpopp/custom-panel/discussions)
- **Documentation**: [Full Documentation Site](https://berntpopp.github.io/custom-panel/)