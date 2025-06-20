# SNP Functionality Documentation

## Overview

The custom-panel tool includes comprehensive SNP (Single Nucleotide Polymorphism) functionality for various genomic applications including sample tracking, ancestry analysis, polygenic risk scoring, and pharmacogenomics. This document covers both the scraping architecture and pipeline usage.

## Table of Contents

1. [SNP Categories](#snp-categories)
2. [Scraping Architecture](#scraping-architecture)
3. [Pipeline Integration](#pipeline-integration)
4. [Configuration](#configuration)
5. [Usage Examples](#usage-examples)
6. [Data Formats](#data-formats)
7. [Troubleshooting](#troubleshooting)

## SNP Categories

The system supports several categories of SNPs:

### 1. Identity SNPs
- **Purpose**: Sample tracking and quality control
- **Use Case**: Verify sample identity across processing steps
- **Sources**: Commercial panels, research publications
- **Examples**: Pengelly panel (24 SNPs), Eurogentest guidelines (14 SNPs)

### 2. Ethnicity/Ancestry SNPs  
- **Purpose**: Population ancestry inference
- **Use Case**: Clinical reporting, population stratification
- **Sources**: Ancestry informative SNPs (AISNPs) databases
- **Examples**: Continental ancestry markers

### 3. Polygenic Risk Score (PRS) SNPs
- **Purpose**: Risk prediction for complex diseases
- **Use Case**: Cancer risk assessment, preventive medicine
- **Sources**: GWAS consortiums, PGS Catalog
- **Examples**: BCAC breast cancer PRS (313 SNPs)

### 4. Deep Intronic ClinVar SNPs
- **Purpose**: Pathogenic variants in non-coding regions
- **Use Case**: Comprehensive variant interpretation
- **Sources**: ClinVar database filtering
- **Examples**: Splice site variants, regulatory variants

### 5. Pharmacogenomics (PGx) SNPs
- **Purpose**: Drug response prediction
- **Use Case**: Personalized medication selection
- **Sources**: PharmGKB, CPIC guidelines
- **Examples**: CYP2C19 variants for clopidogrel response

## Scraping Architecture

### Core Components

#### 1. Base SNP Scraper (`base_snp_scraper.py`)
```python
class BaseSNPScraper:
    """
    Abstract base class for all SNP scrapers.
    Provides common functionality and interface.
    """
    
    def __init__(self, url: str):
        self.url = url
    
    def parse(self) -> Dict[str, Any]:
        """Main parsing method - must be implemented by subclasses"""
        raise NotImplementedError
    
    def validate_rsid(self, rsid: str) -> bool:
        """Validate rsID format (rs followed by digits)"""
        return bool(re.match(r'^rs\d+$', rsid, re.IGNORECASE))
    
    def create_snp_record(self, rsid: str, category: str, panel_specific_name: str) -> Dict:
        """Create standardized SNP record"""
        return {
            "rsid": rsid,
            "source": self.__class__.__name__.replace("Parser", ""),
            "category": category,
            "panel_specific_name": panel_specific_name
        }
```

#### 2. Content Fetching (`downloader.py`)
```python
class PanelDownloader:
    """
    Handles downloading of various file types (PDF, HTML, TSV, etc.)
    with intelligent content type detection and caching.
    """
    
    def download(self, url: str, expected_format: str, method: str = "requests") -> Path:
        """
        Download file from URL with format validation
        
        Args:
            url: Source URL
            expected_format: Expected file type (pdf, html, tsv, etc.)
            method: Download method (requests, selenium)
        
        Returns:
            Path to downloaded file
        """
```

#### 3. Specialized Parsers

Each parser handles specific data sources and formats:

| Parser | Source | Format | SNPs | Technology |
|--------|--------|--------|------|------------|
| `PengellyParser` | Research paper | PDF | 24 | Sample tracking |
| `EurogentestParser` | Guidelines | PDF | 14 | NGS quality control |
| `IDTAmpliconParser` | IDT website | HTML | ~50 | Amplicon sequencing |
| `IDTHybridizationParser` | IDT targets | TSV | 71 | Capture sequencing |
| `AmpliseqParser` | Illumina docs | PDF | 8 | Amplicon sequencing |
| `NimaGenHESTParser` | Product sheet | PDF | 37 | Sample identification |
| `NimaGenHISTParser` | Product sheet | PDF | 34 | Sample tracking |

### Parsing Strategies

#### PDF Parsing
```python
def _extract_rsids_from_pdf(self, pdf_path: Path) -> List[str]:
    """
    Multi-strategy PDF parsing:
    1. Table extraction using pdfplumber
    2. Text extraction with regex patterns  
    3. Fallback to full-text search
    """
    with pdfplumber.open(pdf_path) as pdf:
        # Try table extraction first
        for page in pdf.pages:
            tables = page.extract_tables()
            for table in tables:
                # Extract rsIDs from table cells
                
        # Fallback to text extraction
        for page in pdf.pages:
            text = page.extract_text()
            rsids.extend(re.findall(r'\b(rs\d+)\b', text))
```

#### HTML Parsing
```python
def _extract_rsids_from_html(self, soup: BeautifulSoup) -> List[str]:
    """
    Multi-strategy HTML parsing:
    1. Look for structured tables
    2. Search in specific div classes
    3. Extract from page text as fallback
    """
    # Strategy 1: Structured tables
    for table_class in ['table-condensed', 'table', 'data-table']:
        tables = soup.find_all('table', class_=table_class)
        # Parse table content
        
    # Strategy 2: Specific containers
    for div_class in ['snp-list', 'marker-list']:
        divs = soup.find_all('div', class_=div_class)
        # Extract from div content
```

#### Network Resilience
```python
def parse(self) -> Dict[str, Any]:
    """
    Robust parsing with fallback mechanisms
    """
    try:
        # Primary method: Selenium (for JavaScript content)
        html_content = self.fetch_content(use_selenium=True)
    except Exception as selenium_error:
        logger.warning(f"Selenium failed: {selenium_error}")
        try:
            # Fallback: Standard HTTP requests
            html_content = self.fetch_content(use_selenium=False)
        except Exception as requests_error:
            raise Exception(f"All methods failed: {requests_error}")
```

### Data Standardization

All scrapers produce standardized JSON output:

```json
{
  "panel_name": "pengelly_panel",
  "source_url": "https://example.com/source.pdf",
  "description": "Pengelly et al. (2013) - 56 SNP panel for sample tracking",
  "snps": [
    {
      "rsid": "rs10203363",
      "source": "Pengelly_Panel", 
      "category": "identity",
      "panel_specific_name": "Pengelly et al. 2013"
    }
  ],
  "metadata": {
    "publication": "Genome Medicine 2013",
    "doi": "10.1186/gm492", 
    "source_url": "https://example.com/source.pdf",
    "snp_count": 24
  },
  "retrieval_date": "2025-06-20"
}
```

## Pipeline Integration

### SNP Annotator (`snp_annotator.py`)
```python
class SNPAnnotator:
    """
    Handles SNP processing within the main pipeline
    """
    
    def __init__(self, config: ConfigManager):
        self.config = config
        self.parsers = {
            'simple_rsid': self._parse_simple_rsid,
            'excel_rsid': self._parse_excel_rsid,
            'bcac_prs': self._parse_bcac_prs,
            'manual_csv': self._parse_manual_csv
        }
    
    def process_snps(self, gene_panel: pd.DataFrame) -> Dict[str, Any]:
        """
        Main SNP processing workflow:
        1. Load identity/ethnicity SNPs
        2. Load PRS SNPs  
        3. Extract deep intronic ClinVar variants
        4. Load manual SNPs
        5. Combine and annotate all SNPs
        """
```

### Integration Points

#### 1. After Gene Panel Creation
```python
# In pipeline.py
def run_pipeline():
    # ... gene panel processing ...
    
    if config.snp_processing.enabled:
        snp_annotator = SNPAnnotator(config)
        snp_data = snp_annotator.process_snps(final_gene_panel)
        
        # Add SNP data to outputs
        output_manager.add_snp_data(snp_data)
```

#### 2. Output Generation
```python
# SNP data included in final outputs:
# - snp_panel.xlsx: All SNPs with annotations
# - snp_panel.bed: Genomic coordinates for SNPs
# - snp_summary.html: Interactive SNP browser
```

### Processing Workflow

```mermaid
graph TD
    A[Gene Panel Complete] --> B{SNP Processing Enabled?}
    B -->|Yes| C[Load Identity SNPs]
    B -->|No| Z[Skip SNP Processing]
    
    C --> D[Load Ethnicity SNPs]
    D --> E[Load PRS SNPs]
    E --> F[Extract ClinVar Deep Intronic]
    F --> G[Load Manual SNPs]
    
    G --> H[Combine All SNPs]
    H --> I[Annotate with Coordinates]
    I --> J[Generate SNP Outputs]
    J --> K[Update Final Report]
```

## Configuration

### Main Configuration (`default_config.yml`)

```yaml
snp_processing:
  enabled: true
  
  identity_and_ethnicity:
    enabled: true
    panels:
      - name: "Pengelly_Panel"
        file_path: "data/snp/identity/pengelly_2013.txt"
        parser: "simple_rsid"
      
  prs:
    enabled: true
    panels:
      - name: "BCAC_313_PRS"
        file_path: "data/snp/prs/BCAC_313_PRS.csv"
        parser: "bcac_prs"
        rsid_column: "hm_rsid"
        effect_allele_column: "effect_allele"
        
  deep_intronic_clinvar:
    enabled: true
    vcf_path: "data/clinvar/clinvar.vcf.gz"
    intronic_padding: 50
    
  manual_snps:
    enabled: true
    lists:
      - name: "PGx_CYP2C19"
        file_path: "data/snp/manual/pgx_cyp2c19.csv"
        parser: "manual_csv"
```

### SNP Scrapers Configuration

```yaml
snp_scrapers:
  pengelly_panel:
    enabled: true
    url: "https://genomemedicine.biomedcentral.com/counter/pdf/10.1186/gm492.pdf"
    parser_module: "parse_pengelly"
    parser_class: "PengellyParser"
    output_path: "data/scraped/snp/pengelly_panel.json"
    
  # ... additional scrapers ...
```

## Usage Examples

### 1. Running SNP Scrapers

```bash
# Run all enabled SNP scrapers
cd scrapers
python run_snp_scrapers.py

# Run specific scraper
python run_snp_scrapers.py --scrapers pengelly_panel

# Enable debug logging
python run_snp_scrapers.py --log-level DEBUG
```

### 2. Full Pipeline with SNPs

```bash
# Run complete pipeline including SNP processing
custom-panel run --config config.local.yml

# Run with SNP processing disabled
custom-panel run --config config.yml --disable-snp-processing
```

### 3. SNP-Only Processing

```bash
# Process SNPs without running gene panel pipeline  
custom-panel process-snps --gene-panel results/latest/master_panel.xlsx
```

### 4. Converting Scraped Data

```bash
# Convert scraped JSON to simple text files
custom-panel convert-snp-data \
  --input data/scraped/snp/ \
  --output data/snp/identity/ \
  --format simple_rsid
```

## Data Formats

### Input Formats

#### 1. Simple rsID List
```text
rs10203363
rs1037256
rs10373
rs10883099
```

#### 2. Excel with Metadata
```csv
rs_id,group,population,maf
rs10203363,European,CEU,0.42
rs1037256,African,YRI,0.18
```

#### 3. PRS Format (BCAC)
```csv
hm_rsid,hm_chr,hm_pos,effect_allele,other_allele,effect_weight
rs10069690,5,1279790,T,C,0.0054
rs1011970,9,136154168,T,C,-0.0208
```

#### 4. Manual CSV
```csv
rsID,Gene,Drug,Phenotype,Evidence
rs4244285,CYP2C19,Clopidogrel,Poor metabolizer,1A
rs4986893,CYP2C19,Clopidogrel,Rapid metabolizer,1A
```

### Output Formats

#### 1. Standardized JSON
```json
{
  "panel_name": "combined_identity_snps",
  "snps": [...],
  "metadata": {
    "total_snps": 284,
    "sources": ["Pengelly", "Eurogentest", "IDT", "Illumina"],
    "categories": ["identity", "ethnicity", "prs", "pgx"]
  }
}
```

#### 2. BED Format
```bed
1	1279790	1279791	rs10069690	0	+
9	136154168	136154169	rs1011970	0	+
```

#### 3. Excel Output
- **Sheet 1**: All SNPs with annotations
- **Sheet 2**: Summary by category
- **Sheet 3**: Source breakdown  
- **Sheet 4**: Quality metrics

## Troubleshooting

### Common Issues

#### 1. Network/Proxy Issues
**Problem**: Scrapers fail with proxy errors
```
ERROR: ERR_ACCESS_DENIED from proxy-cbf-2.charite.de
```

**Solution**: System automatically falls back to requests method
```python
# Automatic fallback implemented
try:
    content = fetch_with_selenium(url)
except ProxyError:
    content = fetch_with_requests(url)
```

#### 2. PDF Parsing Issues
**Problem**: Zero SNPs extracted from PDF
```
WARNING: No rsIDs found in expected format
```

**Solution**: Check PDF structure and adjust parsing strategy
```python
# Debug PDF content
with pdfplumber.open(pdf_path) as pdf:
    for page_num, page in enumerate(pdf.pages):
        print(f"Page {page_num}: {page.extract_text()[:200]}")
```

#### 3. Missing Dependencies
**Problem**: Import errors for PDF/HTML parsing
```
ModuleNotFoundError: No module named 'pdfplumber'
```

**Solution**: Install required packages
```bash
pip install pdfplumber beautifulsoup4 selenium webdriver-manager lxml html5lib
```

#### 4. Configuration Issues
**Problem**: SNP processing not running
```
INFO: SNP processing disabled in configuration
```

**Solution**: Enable in configuration
```yaml
snp_processing:
  enabled: true
```

### Debugging Tools

#### 1. Verbose Logging
```bash
python run_snp_scrapers.py --log-level DEBUG
```

#### 2. Single Scraper Testing
```python
from scrapers.snp_parsers.parse_pengelly import PengellyParser

parser = PengellyParser("https://example.com/paper.pdf")
result = parser.parse()
print(f"Found {len(result['snps'])} SNPs")
```

#### 3. Configuration Validation
```bash
custom-panel validate-config --config config.local.yml
```

#### 4. Cache Management
```bash
# Clear download cache
rm -rf data/snp/downloads/

# Clear API cache  
rm -rf .cache/
```

### Performance Optimization

#### 1. Concurrent Processing
```python
# Process multiple scrapers in parallel
with ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(scraper.parse) for scraper in scrapers]
```

#### 2. Caching Strategy
```python
# Cache downloaded files
@lru_cache(maxsize=100)
def download_and_cache(url: str) -> Path:
    return downloader.download(url)
```

#### 3. Batch Processing
```python
# Process SNPs in batches
for batch in batch_iterator(snps, batch_size=1000):
    annotated_batch = annotate_snps(batch)
```

## Advanced Usage

### Custom Parser Development

#### 1. Create New Parser
```python
class CustomSNPParser(BaseSNPScraper):
    def parse(self) -> Dict[str, Any]:
        # Download source data
        data = self.fetch_content()
        
        # Extract rsIDs using custom logic
        rsids = self._extract_rsids(data)
        
        # Create SNP records
        snps = [self.create_snp_record(rsid, "custom", "Custom Panel") 
                for rsid in rsids]
        
        return {
            "panel_name": "custom_panel",
            "source_url": self.url,
            "snps": snps,
            "metadata": {"snp_count": len(snps)}
        }
```

#### 2. Register Parser
```yaml
# Add to config.yml
snp_scrapers:
  custom_panel:
    enabled: true
    url: "https://example.com/data"
    parser_module: "parse_custom"
    parser_class: "CustomSNPParser"
    output_path: "data/scraped/snp/custom_panel.json"
```

### Integration with External Tools

#### 1. PLINK Integration
```bash
# Generate PLINK files from SNP panel
custom-panel export-plink \
  --snp-panel results/snp_panel.xlsx \
  --output results/snp_panel.plink
```

#### 2. VCF Generation
```bash
# Create VCF template for SNP panel
custom-panel export-vcf \
  --snp-panel results/snp_panel.xlsx \
  --reference GRCh38 \
  --output results/snp_panel.vcf
```

#### 3. Database Integration
```python
# Load SNP data into database
from custom_panel.database import SNPDatabase

db = SNPDatabase("postgresql://user:pass@host/db")
db.load_snp_panel("results/snp_panel.xlsx")
```

This comprehensive documentation covers all aspects of the SNP functionality, from the scraping architecture to pipeline integration and practical usage examples.