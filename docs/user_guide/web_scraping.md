# Web Scraping System

Custom Panel includes a comprehensive web scraping framework for extracting gene lists from commercial diagnostic panel websites. This replaces manual data collection with an automated, maintainable system.

## Overview

The scraping system:
- **Automates** collection from 14 commercial providers
- **Standardizes** output to consistent JSON format
- **Validates** gene symbols with quality checks
- **Isolates** failures to prevent cascade effects
- **Maintains** compatibility with changing websites

## Architecture

```
scrapers/
├── run_scrapers.py          # Master runner with CLI
├── parsers/                 # Individual parser implementations
│   ├── base_parser.py       # Abstract base class
│   ├── parse_myriad.py      # Myriad Genetics parser
│   ├── parse_blueprint.py   # Blueprint Genetics parser
│   └── ...                  # 14 total parsers
└── README.md               # Scraper-specific documentation
```

### Design Principles

1. **Decoupled**: Scrapers run independently from main tool
2. **Fault Tolerant**: Individual failures don't affect others  
3. **Maintainable**: Easy to update for website changes
4. **Consistent**: Standardized output format
5. **Extensible**: Simple to add new providers

## Running the Scrapers

### Standalone Execution

```bash
# Run all enabled scrapers
python scrapers/run_scrapers.py

# Run specific scrapers only
python scrapers/run_scrapers.py --names myriad_myrisk blueprint_genetics

# Preview what would be executed
python scrapers/run_scrapers.py --dry-run

# Custom output directory
python scrapers/run_scrapers.py --output-dir /path/to/output

# Enable verbose logging
python scrapers/run_scrapers.py --verbose
```

### Integration with Custom Panel

```bash
# Fetch pre-scraped commercial panel data
custom-panel fetch commercial_panels --output-dir results
```

## Currently Implemented Scrapers

### 1. Myriad Genetics (`parse_myriad.py`)
- **URL**: https://myriad.com/gene-table/
- **Panel**: myRisk Hereditary Cancer Panel
- **Method**: Static HTML with BeautifulSoup
- **Genes**: ~35 high-risk cancer genes

### 2. Blueprint Genetics (`parse_blueprint.py`)
- **URL**: Multiple sub-panels (19 total)
- **Panels**: Comprehensive hereditary cancer panels
- **Method**: Dynamic content with Selenium
- **Genes**: ~200 genes across all panels

### 3. Invitae (`parse_invitae.py`)
- **URL**: https://www.invitae.com/en/providers/test-catalog/test-01101
- **Panel**: Multi-Cancer Panel
- **Method**: JavaScript-rendered with Selenium
- **Genes**: ~80 cancer predisposition genes

### 4. GeneDx (`parse_genedx.py`)
- **URL**: https://www.genedx.com/tests/detail/oncogenedx-custom-panel-871
- **Panel**: Comprehensive Cancer Panel
- **Method**: Static HTML parsing
- **Genes**: ~70 hereditary cancer genes

### 5. Fulgent Genetics (`parse_fulgent.py`)
- **URL**: https://www.fulgentgenetics.com/comprehensivecancer-full
- **Panel**: Comprehensive Cancer Panel
- **Method**: BeautifulSoup with fallback strategies
- **Genes**: ~130 cancer genes

### Additional Providers

6. **Centogene** - Solid tumor panel
7. **CEGAT** - Tumor syndrome panel
8. **MGZ Munich** - German cancer panel
9. **University of Chicago** - Academic cancer panel
10. **Prevention Genetics** - Hereditary cancer panel
11. **ARUP Laboratories** - Clinical cancer panel
12. **Cincinnati Children's** - Pediatric focus
13. **NeoGenomics** - Oncology specialists
14. **Natera** - Hereditary cancer test

## Parser Implementation

### Base Parser Class

All parsers inherit from `BaseParser`:

```python
class BaseParser(ABC):
    """Abstract base parser for commercial panels."""
    
    def __init__(self, config: dict[str, Any]):
        self.name = config["name"]
        self.url = config["url"]
        self.output_path = config["output_path"]
        
    @abstractmethod
    def parse(self) -> list[str]:
        """Extract gene symbols from website."""
        pass
        
    def save_results(self, genes: list[str]) -> None:
        """Save standardized JSON output."""
        # Implemented in base class
```

### Parsing Strategies

#### Static HTML (BeautifulSoup)

For simple HTML pages:

```python
def parse(self) -> list[str]:
    response = requests.get(self.url, timeout=30)
    soup = BeautifulSoup(response.content, "html.parser")
    
    # Primary: XPath-like selection
    genes = soup.select("td.gene-symbol")
    
    # Fallback: Text pattern matching
    if not genes:
        genes = soup.find_all(text=re.compile(r'^[A-Z][A-Z0-9]+$'))
        
    return self._clean_gene_list(genes)
```

#### Dynamic JavaScript (Selenium)

For JavaScript-rendered content:

```python
def parse(self) -> list[str]:
    options = webdriver.ChromeOptions()
    options.add_argument("--headless")
    
    with webdriver.Chrome(options=options) as driver:
        driver.get(self.url)
        
        # Wait for content to load
        wait = WebDriverWait(driver, 10)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "gene-list")))
        
        # Extract genes
        elements = driver.find_elements(By.CSS_SELECTOR, ".gene-name")
        return [elem.text.strip() for elem in elements]
```

### Gene Validation

All parsers apply consistent validation:

```python
def _clean_and_validate(self, genes: list[str]) -> list[str]:
    """Clean and validate gene symbols."""
    cleaned = []
    
    for gene in genes:
        # Remove common suffixes
        gene = re.sub(r'\*|\(.*\)|\s+', '', gene)
        
        # Validate format
        if self._is_valid_gene(gene):
            cleaned.append(gene)
            
    return list(dict.fromkeys(cleaned))  # Remove duplicates

def _is_valid_gene(self, symbol: str) -> bool:
    """Check if string is likely a gene symbol."""
    return (
        1 <= len(symbol) <= 20 and
        re.match(r'^[A-Z][A-Z0-9\-\_]*$', symbol) and
        symbol not in SKIP_TERMS
    )
```

## Output Format

All scrapers produce standardized JSON:

```json
{
  "panel_name": "myriad_myrisk",
  "source_url": "https://myriad.com/gene-table/",
  "retrieval_date": "2024-01-15",
  "gene_count": 35,
  "genes": [
    "APC",
    "ATM", 
    "BRCA1",
    "BRCA2",
    "CDH1",
    "CHEK2",
    "MLH1",
    "MSH2",
    "MSH6",
    "PALB2",
    "PMS2",
    "PTEN",
    "STK11",
    "TP53"
  ]
}
```

## Adding New Scrapers

### Step 1: Create Parser Class

Create `scrapers/parsers/parse_newprovider.py`:

```python
from typing import Any
from .base_parser import BaseParser

class NewProviderParser(BaseParser):
    """Parser for NewProvider gene panel."""
    
    def parse(self) -> list[str]:
        """Extract genes from NewProvider website."""
        # Implement scraping logic
        response = self._fetch_page()
        genes = self._extract_genes(response)
        return self._clean_gene_list(genes)
```

### Step 2: Add Configuration

Update `custom_panel/config/default_config.yml`:

```yaml
scrapers:
  newprovider:
    enabled: true
    url: "https://newprovider.com/panel"
    parser_module: "parse_newprovider"
    parser_class: "NewProviderParser"
    output_path: "data/scraped/newprovider.json"
```

### Step 3: Add to Commercial Panels

```yaml
data_sources:
  Commercial_Panels:
    panels:
      - name: "NewProvider"
        file_path: "data/scraped/newprovider.json"
        evidence_score: 0.8
```

### Step 4: Test the Parser

```bash
# Test individual parser
python scrapers/run_scrapers.py --names newprovider --verbose

# Verify output
cat data/scraped/newprovider.json
```

## Error Handling

The system includes robust error handling:

### Network Errors
- Timeout after 30 seconds
- Retry logic with exponential backoff
- Graceful failure with error logging

### Parsing Errors
- Multiple fallback strategies
- Pattern-based extraction as last resort
- Detailed error messages for debugging

### Validation Errors
- Skip invalid gene symbols
- Log suspicious patterns
- Continue with valid genes

## Maintenance

### Monitoring Website Changes

1. **Regular Testing**: Run scrapers weekly to detect changes
2. **Gene Count Validation**: Alert if count changes significantly
3. **Visual Inspection**: Periodically verify scraped genes
4. **Error Tracking**: Monitor logs for parsing failures

### Updating Parsers

When websites change:

1. **Identify Changes**: Compare HTML structure
2. **Update Selectors**: Modify CSS/XPath selectors
3. **Add Fallbacks**: Implement alternative strategies
4. **Test Thoroughly**: Verify gene extraction
5. **Document Changes**: Note modifications in code

### Common Issues

**JavaScript Loading**
- Increase wait times
- Add explicit waits for elements
- Use network idle conditions

**Authentication Required**
- Add login automation
- Use session cookies
- Consider API alternatives

**Rate Limiting**
- Add delays between requests
- Rotate user agents
- Respect robots.txt

## Best Practices

1. **Respect Websites**: Follow robots.txt and terms of service
2. **Cache Results**: Don't re-scrape unnecessarily  
3. **Version Control**: Track scraper changes carefully
4. **Test Regularly**: Automated tests for each parser
5. **Document Thoroughly**: Clear comments for maintainers

## Integration with Pipeline

The scraped data integrates seamlessly:

```bash
# 1. Run scrapers (periodic update)
python scrapers/run_scrapers.py

# 2. Use in pipeline
custom-panel run --output-dir results

# 3. View commercial panel contribution
grep "Commercial_Panels" results/run_*/run_summary.json
```

## Next Steps

- [Data Sources Overview](./data_sources.md) - Understanding all data sources
- [Configuration Guide](./configuration.md) - Customize commercial panel settings
- [Scoring System](./scoring_system.md) - How commercial panels affect scores