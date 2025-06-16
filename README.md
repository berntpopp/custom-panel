# Custom Panel

A modern Python tool for gene panel curation and aggregation from multiple genomic databases.

[![CI](https://github.com/your-username/custom-panel/workflows/CI/badge.svg)](https://github.com/your-username/custom-panel/actions)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Custom Panel is a comprehensive tool for creating, managing, and curating gene panels for clinical genomics applications. It aggregates gene information from multiple trusted sources, applies configurable scoring algorithms, and generates standardized output formats suitable for downstream analysis.

### Key Features

- **Multi-source data aggregation**: Integrates data from PanelApp, ACMG recommendations, in-house panels, manual curation lists, and HPO/OMIM
- **Manual curation support**: Process custom gene lists from Excel, CSV, and text files with configurable parameters
- **HPO/OMIM integration**: Automatic identification of neoplasm-associated genes using HPO ontology and OMIM data
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

# Fetch manual curation lists
custom-panel fetch manual --output-dir results/manual

# Fetch HPO neoplasm genes
custom-panel fetch hpo --output-dir results/hpo
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
      manual_curation: 1.1
      hpo_neoplasm: 0.7
    somatic:
      cosmic: 1.2
      commercial_panels: 0.9
      manual_curation: 1.0
      hpo_neoplasm: 0.9
  
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
  SOURCE  Data source (panelapp, inhouse, acmg, manual, hpo)

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
4. **Manual Curation**: Custom gene lists from literature review and expert curation
5. **HPO/OMIM Neoplasm**: Automatic identification of cancer-associated genes using HPO ontology
6. **COSMIC**: Cancer Gene Census (planned)
7. **Commercial Panels**: Web scraping of 15 commercial diagnostic panel websites

### Commercial Panel Web Scrapers

The tool includes a comprehensive web scraping framework for extracting gene lists from commercial diagnostic panel websites. This replaces brittle manual processes with an automated, maintainable system.

#### Supported Commercial Panels

The scrapers framework currently supports **15 commercial diagnostic panel providers**:

1. **Myriad Genetics** - myRisk Hereditary Cancer Panel
2. **Blueprint Genetics** - Hereditary Cancer Panels (19 sub-panels)
3. **Centogene** - Solid Tumor Panel
4. **Fulgent Genetics** - Comprehensive Cancer Panel
5. **CEGAT** - Tumor Syndromes Panel
6. **Invitae** - Comprehensive Cancer Panel
7. **MGZ Munich** - Hereditary Cancer Panel
8. **University of Chicago** - Hereditary Cancer Panel
9. **Prevention Genetics** - Cancer Panel
10. **Mayo Clinic Labs** - Hereditary Expanded Cancer Panel
11. **GeneDx** - Comprehensive Cancer Panel
12. **ARUP Laboratories** - Hereditary Cancer Panel
13. **Cincinnati Children's (CGL)** - Hereditary Cancer Panel
14. **NeoGenomics** - Comprehensive Cancer Panel
15. **Natera** - Hereditary Cancer Test

#### Running the Scrapers

**Standalone Scraper Runner:**
```bash
# Run all enabled scrapers
python scrapers/run_scrapers.py

# Run specific scrapers
python scrapers/run_scrapers.py --names myriad_myrisk blueprint_genetics

# Dry run to see what would be executed
python scrapers/run_scrapers.py --dry-run

# Run with custom output directory
python scrapers/run_scrapers.py --output-dir /path/to/custom/output

# Enable verbose logging
python scrapers/run_scrapers.py --verbose
```

**Integrated with Custom Panel CLI:**
```bash
# Fetch commercial panel data (requires scraped JSON files)
custom-panel fetch commercial_panels --output-dir results/commercial
```

#### Architecture Overview

**Decoupled Design:**
- **Independent Scrapers**: Run separately from the main custom-panel tool
- **Standardized Output**: All scrapers produce consistent JSON format
- **Fault Isolation**: Scraping failures don't affect core tool functionality
- **Easy Maintenance**: Update scrapers without touching core logic

**Framework Components:**
```
scrapers/
├── run_scrapers.py          # Master runner script with CLI
├── parsers/                 # Individual parser implementations
│   ├── base_parser.py       # Abstract base class
│   ├── parse_myriad.py      # Myriad Genetics parser
│   ├── parse_blueprint.py   # Blueprint Genetics parser
│   └── ...                  # 15 total parsers
└── README.md               # Scraper-specific documentation
```

#### Parser Implementation Details

**Two Parsing Approaches:**
1. **Static Content** (`requests` + `BeautifulSoup`): For simple HTML pages
2. **Dynamic Content** (`Selenium` + `ChromeDriver`): For JavaScript-heavy sites

**Parsing Strategies:**
- **XPath Selectors**: Precise element targeting (matches original R implementation)
- **CSS Selectors**: Modern web scraping approach
- **Regex Extraction**: Fallback for complex content patterns
- **Multi-level Fallbacks**: Robust error handling with alternative strategies

**Example Parser Structure:**
```python
class MyriadParser(BaseParser):
    def parse(self) -> list[str]:
        # Primary: XPath selector matching R script
        genes = self._parse_with_xpath('//td[contains(@class,"gene")]')
        
        # Fallback: CSS selector approach  
        if not genes:
            genes = self._parse_with_css('td[class*="gene"]')
            
        # Final fallback: Regex extraction
        if not genes:
            genes = self._parse_with_regex(r'\b[A-Z][A-Z0-9]{1,7}\b')
            
        return self._clean_and_validate(genes)
```

#### Output Format

**Standardized JSON Schema:**
```json
{
  "panel_name": "myriad_myrisk",
  "source_url": "https://myriad.com/gene-table/",
  "retrieval_date": "2024-01-15",
  "genes": ["BRCA1", "BRCA2", "TP53", "PTEN", "ATM", ...]
}
```

**Gene Processing:**
- **Symbol Cleaning**: Remove asterisks, parenthetical content, whitespace
- **Validation**: Length (1-20 chars), character set (A-Z, 0-9, -, _), contains letters
- **Skip Terms**: Filters out common non-genes ("GENE", "DNA", "PANEL", etc.)
- **Deduplication**: Removes duplicates while preserving order
- **HGNC Integration**: Validates against HGNC gene nomenclature (in main tool)

#### Configuration

**Scraper Configuration** (`custom_panel/config/default_config.yml`):
```yaml
scrapers:
  myriad_myrisk:
    enabled: true
    url: "https://myriad.com/gene-table/"
    parser_module: "parse_myriad"
    parser_class: "MyriadParser"
    output_path: "data/scraped/myriad_myrisk.json"
  
  blueprint_genetics:
    enabled: true
    url: "https://blueprintgenetics.com/tests/panels/hereditary-cancer/"
    subpanel_urls:
      - "https://blueprintgenetics.com/tests/panels/hematology/comprehensive-hematology-and-hereditary-cancer-panel/"
      - "https://blueprintgenetics.com/tests/panels/hereditary-cancer/comprehensive-hereditary-cancer-panel/"
      # ... 17 more sub-panel URLs
    parser_module: "parse_blueprint"
    parser_class: "BlueprintParser"
    output_path: "data/scraped/blueprint_genetics.json"
```

**Integration Configuration:**
```yaml
data_sources:
  commercial_panels:
    enabled: true
    panels:
      - name: "Myriad_myRisk"
        file_path: "data/scraped/myriad_myrisk.json"
        evidence_score: 0.8
      - name: "Blueprint_Genetics"
        file_path: "data/scraped/blueprint_genetics.json"
        evidence_score: 0.8
```

#### Maintenance and Updates

**Website Change Monitoring:**
- Scrapers include robust fallback strategies for layout changes
- XPath selectors aligned with original R implementation for consistency
- CSS selectors provide modern alternative approaches
- Comprehensive logging for debugging failed extractions

**Adding New Providers:**
1. Create new parser class inheriting from `BaseParser`
2. Implement `parse()` method with site-specific logic
3. Add configuration entry to `default_config.yml`
4. Add integration to commercial panels data source
5. Include comprehensive tests with mock HTML fixtures

**Testing Framework:**
```bash
# Run scraper tests
poetry run pytest tests/test_scrapers.py -v

# Test specific parser
poetry run pytest tests/test_scrapers.py::TestParsers::test_myriad_parser

# Test CLI integration
poetry run pytest tests/test_cli_integration.py -v
```

#### Error Handling and Monitoring

**Robust Error Handling:**
- Network timeouts and retry logic
- Graceful degradation with multiple parsing strategies
- Comprehensive logging for debugging
- Selenium driver management and cleanup

**Monitoring:**
- Success/failure metrics per scraper
- Gene count validation against expected ranges
- Output file validation and formatting checks
- Performance metrics and execution timing

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

### Manual Curation Lists

Process manually curated gene lists with flexible configuration:

```yaml
data_sources:
  manual_curation:
    enabled: true
    lists:
      - name: "Literature_Review_2024"
        file_path: "data/manual/literature_genes_2024.txt"
        gene_column: "gene_symbol"  # Column containing gene symbols
        evidence_score: 1.2
      - name: "Expert_Panel_Recommendations"
        file_path: "data/manual/expert_panel_genes.csv"
        gene_column: "Gene"
        evidence_score: 1.1
      - name: "Clinical_Variant_Database"
        file_path: "data/manual/clinical_variants.xlsx"
        gene_column: "Gene Symbol"
        evidence_score: 1.0
```

Features:
- Support for multiple file formats (Excel, CSV, TXT)
- Configurable gene column names
- Individual evidence scoring per list
- Automatic gene extraction and standardization

### HPO/OMIM Neoplasm Integration

Automatically identify cancer-associated genes using HPO ontology:

```yaml
data_sources:
  hpo_neoplasm:
    enabled: true
    # Automatically search for neoplasm-related HPO terms
    use_neoplasm_search: true
    # Specific HPO terms to include
    specific_hpo_terms:
      - "HP:0002664"  # Neoplasm
      - "HP:0030358"  # Non-small cell lung carcinoma
      - "HP:0002860"  # Squamous cell carcinoma
    # Optional OMIM file for additional associations
    omim_file_path: "data/hpo/omim_neoplasm_genes.csv"
    omim_gene_column: "Gene Symbol"
    omim_disease_column: "Disease"
    omim_id_column: "OMIM ID"
    # Scoring parameters
    base_evidence_score: 0.6
    max_hierarchy_depth: 5
```

Features:
- Automatic neoplasm term discovery in HPO
- Hierarchical HPO term traversal
- OMIM file integration for additional gene-disease associations
- Configurable evidence scoring based on data source diversity
- Comprehensive gene-phenotype relationship mapping

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
- **HPO** (Human Phenotype Ontology) for phenotype-gene associations
- **OMIM** for comprehensive gene-disease relationship data

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
