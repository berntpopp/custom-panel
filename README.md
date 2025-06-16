# Custom Panel

A modern Python tool for gene panel curation and aggregation from multiple genomic databases.

[![CI](https://github.com/your-username/custom-panel/workflows/CI/badge.svg)](https://github.com/your-username/custom-panel/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Custom Panel is a comprehensive tool for creating, managing, and curating gene panels for clinical genomics applications. It aggregates gene information from multiple trusted sources, applies configurable scoring algorithms, and generates standardized output formats suitable for downstream analysis.

### Key Features

- **Multi-source data aggregation**: Integrates data from PanelApp, ACMG recommendations, in-house panels, and more
- **Intelligent scoring system**: Configurable evidence weighting and decision thresholds
- **Gene standardization**: Automatic gene symbol standardization using HGNC
- **Genomic annotation**: Rich annotation with Ensembl coordinates, transcripts, and MANE information
- **Flexible output formats**: Excel, CSV, Parquet, and BED file generation
- **Modern architecture**: Built with Python 3.10+, Poetry, and type hints
- **Comprehensive testing**: Full test suite with pytest
- **CI/CD ready**: GitHub Actions workflows for testing and linting

## Installation

### Using Poetry (Recommended)

```bash
# Clone the repository
git clone https://github.com/your-username/custom-panel.git
cd custom-panel

# Install with Poetry
poetry install

# Activate the virtual environment
poetry shell
```

### Using pip

```bash
# Clone the repository
git clone https://github.com/your-username/custom-panel.git
cd custom-panel

# Install in development mode
pip install -e .
```

### From PyPI (when published)

```bash
pip install custom-panel
```

## Quick Start

### 1. Check Configuration

```bash
custom-panel config-check
```

### 2. Run the Complete Pipeline

```bash
custom-panel run --output-dir results
```

### 3. Fetch Data from Specific Sources

```bash
# Fetch from PanelApp
custom-panel fetch panelapp --output-dir results/panelapp

# Fetch in-house panels
custom-panel fetch inhouse --output-dir results/inhouse

# Fetch ACMG incidental findings
custom-panel fetch acmg --output-dir results/acmg
```

### 4. Search Available Panels

```bash
custom-panel search-panels "cancer"
```

## Configuration

The tool uses a YAML configuration file to define data sources, scoring parameters, and output options. A comprehensive default configuration is provided at `custom_panel/config/default_config.yml`.

### Key Configuration Sections

#### Data Sources

```yaml
data_sources:
  panelapp:
    enabled: true
    panels:
      - name: "PanelApp_UK"
        base_url: "https://panelapp.genomicsengland.co.uk/api/v1/panels"
        panels:
          - id: 1
            name: "Childhood_solid_tumours"
    evidence_scores:
      "Green": 1.0
      "Amber": 0.5
      "Red": 0.1
  
  inhouse_panels:
    enabled: true
    panels:
      - name: "Leipzig-Cancer_v6"
        file_path: "data/Leipzig-Cancer_v6_gene-list.xlsx"
        gene_column: "Gene_Symbol"
        evidence_score: 1.0
```

#### Scoring Configuration

```yaml
scoring:
  category_weights:
    germline:
      panelapp: 1.0
      acmg_incidental: 1.5
      inhouse_panels: 1.2
    somatic:
      cosmic: 1.2
      commercial_panels: 0.9
  
  thresholds:
    germline_threshold: 2.0
    somatic_threshold: 1.5
    min_sources: 1
```

## Command Line Interface

### Main Commands

#### `run` - Execute the complete pipeline

```bash
custom-panel run [OPTIONS]

Options:
  -c, --config-file TEXT        Configuration file path
  -o, --output-dir TEXT         Output directory
  --germline-threshold FLOAT    Override germline score threshold
  --somatic-threshold FLOAT     Override somatic score threshold
  --log-level TEXT              Log level (DEBUG, INFO, WARNING, ERROR)
  --dry-run                     Run without generating output files
```

#### `fetch` - Fetch data from a single source

```bash
custom-panel fetch SOURCE [OPTIONS]

Arguments:
  SOURCE  Data source (panelapp, inhouse, acmg)

Options:
  -c, --config-file TEXT  Configuration file path
  -o, --output-dir TEXT   Output directory
  -f, --format TEXT       Output format (parquet, csv, excel)
```

#### `config-check` - Validate configuration

```bash
custom-panel config-check [OPTIONS]

Options:
  -c, --config-file TEXT  Configuration file path
```

#### `search-panels` - Search PanelApp panels

```bash
custom-panel search-panels QUERY [OPTIONS]

Arguments:
  QUERY  Search term for panel names

Options:
  -c, --config-file TEXT  Configuration file path
```

## Data Sources

### Supported Sources

1. **PanelApp**: UK Genomics England and Australia PanelApp APIs
2. **In-house Panels**: Local Excel, CSV, or text files
3. **ACMG Incidental Findings**: ACMG SF v3.2 recommendations
4. **COSMIC**: Cancer Gene Census (planned)
5. **Commercial Panels**: PDF and Excel parsing (planned)
6. **HPO**: Human Phenotype Ontology associations (planned)

### Adding Custom Data Sources

To add your own in-house panels, update the configuration:

```yaml
data_sources:
  inhouse_panels:
    enabled: true
    panels:
      - name: "My_Custom_Panel"
        file_path: "data/my_genes.xlsx"
        gene_column: "Gene_Symbol"
        sheet_name: "Sheet1"  # Optional for Excel files
        evidence_score: 1.0
```

Supported file formats:
- Excel (.xlsx, .xls)
- CSV (.csv)
- TSV (.tsv)
- Plain text (.txt) - one gene per line

## Output Formats

### Generated Files

- **master_panel.xlsx**: Comprehensive Excel file with all data
- **master_panel.csv**: CSV format for programmatic use
- **master_panel.parquet**: Efficient binary format
- **germline_panel.bed**: BED file for germline genes
- **somatic_panel.bed**: BED file for somatic genes
- **combined_panel.bed**: BED file for all included genes

### Data Schema

The standardized output includes:

| Column | Description |
|--------|-------------|
| `approved_symbol` | HGNC-approved gene symbol |
| `hgnc_id` | HGNC identifier |
| `germline_score` | Weighted germline evidence score |
| `somatic_score` | Weighted somatic evidence score |
| `total_score` | Combined evidence score |
| `source_count` | Number of supporting sources |
| `source_names` | Semicolon-separated source names |
| `include_germline` | Boolean inclusion for germline panels |
| `include_somatic` | Boolean inclusion for somatic panels |
| `chromosome` | Chromosome location |
| `gene_start` | Gene start position (1-based) |
| `gene_end` | Gene end position |
| `gene_size` | Gene size in base pairs |
| `canonical_transcript` | Ensembl canonical transcript ID |
| `mane_transcript` | MANE transcript ID |

## Development

### Setting up Development Environment

```bash
# Clone and install in development mode
git clone https://github.com/your-username/custom-panel.git
cd custom-panel
poetry install

# Install pre-commit hooks
poetry run pre-commit install

# Run tests
poetry run pytest

# Run all linting and formatting
./scripts/lint.sh

# Or run individual tools
poetry run ruff check .           # Linting
poetry run ruff format .          # Formatting  
poetry run ruff check . --select I001  # Import sorting
poetry run mypy custom_panel/     # Type checking
```

### Running Tests

```bash
# Run all tests
poetry run pytest

# Run with coverage
poetry run pytest --cov=custom_panel

# Run specific test file
poetry run pytest tests/test_core.py

# Run with verbose output
poetry run pytest -v
```

### Code Style

This project uses:
- **Ruff** for linting and formatting
- **MyPy** for type checking
- **Black** code style (via Ruff)
- **isort** import sorting (via Ruff)

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Run the test suite (`poetry run pytest`)
6. Run linting and formatting (`./scripts/lint.sh`)
7. Commit your changes (`git commit -m 'Add amazing feature'`)
8. Push to the branch (`git push origin feature/amazing-feature`)
9. Open a Pull Request

## Scientific Background

### Methodology

The tool implements a weighted evidence aggregation system:

1. **Source Integration**: Data from multiple sources is standardized to a common schema
2. **Gene Standardization**: All gene symbols are validated and standardized using HGNC
3. **Evidence Scoring**: Each source contributes evidence scores weighted by reliability and category
4. **Decision Logic**: Configurable thresholds determine gene inclusion in final panels
5. **Genomic Annotation**: Rich annotation with current genomic coordinates and transcript information

### Quality Control

- HGNC validation ensures gene symbol accuracy
- Duplicate removal based on standardized symbols
- Coordinate validation through Ensembl
- Configurable minimum evidence thresholds
- Source provenance tracking

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **HGNC** for gene nomenclature standardization
- **Ensembl** for genomic coordinates and annotation
- **PanelApp** (UK Genomics England) for curated gene panels
- **ACMG** for incidental findings recommendations

## Citation

If you use this tool in your research, please cite:

```
Custom Panel: A modern Python tool for gene panel curation and aggregation
[Your Name et al.]
[Journal/Repository] [Year]
```

## Support

- **Issues**: [GitHub Issues](https://github.com/your-username/custom-panel/issues)
- **Discussions**: [GitHub Discussions](https://github.com/your-username/custom-panel/discussions)
- **Documentation**: [Wiki](https://github.com/your-username/custom-panel/wiki)
