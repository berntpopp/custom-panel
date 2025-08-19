# Custom Panel

[![CI](https://github.com/berntpopp/custom-panel/actions/workflows/ci.yml/badge.svg)](https://github.com/berntpopp/custom-panel/actions/workflows/ci.yml)
[![Documentation](https://github.com/berntpopp/custom-panel/actions/workflows/docs.yml/badge.svg)](https://github.com/berntpopp/custom-panel/actions/workflows/docs.yml)
[![codecov](https://codecov.io/gh/berntpopp/custom-panel/branch/main/graph/badge.svg)](https://codecov.io/gh/berntpopp/custom-panel)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://github.com/berntpopp/custom-panel)
[![Code style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Type checked: mypy](https://img.shields.io/badge/type%20checked-mypy-blue.svg)](https://github.com/python/mypy)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Docs](https://img.shields.io/badge/docs-view%20site-blue.svg)](https://berntpopp.github.io/custom-panel/)

A modern Python tool for gene panel curation and aggregation from multiple genomic databases.

`custom-panel` is a comprehensive bioinformatics tool for creating, managing, and curating gene panels. It aggregates gene information from trusted sources, applies a configurable scoring system, and generates standardized, analysis-ready output formats.

## Key Features

- **Multi-source data aggregation**: Integrates data from PanelApp, ACMG recommendations, in-house panels, manual curation lists, ClinGen, TheGenCC, COSMIC, and HPO/OMIM
- **Intelligent scoring system**: Configurable evidence weighting with veto capabilities for critical sources
- **Gene standardization**: Automatic gene symbol standardization using HGNC
- **Genomic annotation**: Rich annotation with Ensembl coordinates, transcripts, and MANE information
- **Flexible output formats**: Excel, CSV, Parquet, and BED file generation
- **Modern architecture**: Built with Python 3.10+, uv, and comprehensive type hints

## Quick Start

```bash
# Set up development environment (installs all dependencies)
make dev-setup

# Run the complete pipeline
make run ARGS="run --output-dir results"

# Check configuration
make config-check

# Development workflow
make help       # Show all available commands
make quality    # Run linting and type checking
make test       # Run tests
make docs-serve # Serve documentation locally
```

### Installation Options

```bash
# Install main dependencies only
make install

# Install all dependencies including dev tools
make install-dev

# Install with scrapers dependencies (for web scraping)
make install-scrapers
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