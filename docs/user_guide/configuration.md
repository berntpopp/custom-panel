# Configuration Guide

Learn how to configure Custom Panel's scoring system, data sources, and output options.

## Configuration Overview

Custom Panel uses YAML configuration files to control:
- Data source selection and weighting
- Gene scoring algorithm parameters
- Output formats and quality control
- Performance and caching settings

## Configuration Files

### Default Configuration

The main configuration is in `custom_panel/config/default_config.yml`. This contains all available options with sensible defaults.

### Local Configuration

Create `config.local.yml` to override defaults without modifying the main config:

```yaml
# config.local.yml - overrides for local environment
data_sources:
  COSMIC_Germline:
    # Add your COSMIC credentials here
    username: "your_email@example.com"
    password: "your_password"

  HPO_Neoplasm:
    # Add your OMIM API key
    omim_genemap2_url: "https://data.omim.org/downloads/YOUR_TOKEN/genemap2.txt"
```

## Gene Scoring System

Custom Panel implements a comprehensive scoring system that balances clinical utility with evidence quality.

### Core Scoring Formula

```
Final Score = Σ(source_evidence_score × internal_confidence_score × source_group_weight)
```

### Source Evidence Scores (0.0-1.5)

Evidence scores reflect the clinical reliability of each data source:

```yaml
data_sources:
  ACMG_Incidental_Findings:
    evidence_score: 1.5  # Highest - evidence-based clinical guidelines
  
  ClinGen:
    evidence_score: 1.2  # High quality gene-disease validity
    
  Inhouse_Panels:
    evidence_score: 1.2  # High trust for clinical laboratory use
    
  Manual_Curation:
    evidence_score: 1.2  # Expert reviewed
    
  PanelApp:
    evidence_score: 1.0  # Community consensus standard
    
  COSMIC_Germline:
    evidence_score: 0.9  # Established cancer gene resource
    
  Commercial_Panels:
    evidence_score: 0.8  # Variable quality control
    
  HPO_Neoplasm:
    evidence_score: 0.7  # Automated associations
```

### Classification Multipliers

High-quality sources include evidence classification systems:

**ClinGen/TheGenCC Classification Scores:**
```yaml
classification_scores:
  "Definitive": 1.5
  "Strong": 1.2
  "Moderate": 1.0
  "Limited": 0.3
  "Disputed": 0.3
  "Refuted": 0.0  # Excluded
```

**PanelApp Evidence Levels:**
```yaml
evidence_scores:
  "3": 1.0  # Green - high confidence
  "2": 0.5  # Amber - moderate confidence
  "1": 0.0  # Red - explicitly not recommended
```

### Source Group Weights (Final Multipliers)

Final weighting reflects source priority in clinical decision-making:

```yaml
scoring:
  source_group_weights:
    ACMG_Incidental_Findings: 1.5  # Highest priority
    ClinGen: 1.2
    TheGenCC: 1.2
    Inhouse_Panels: 1.2
    Manual_Curation: 1.1
    PanelApp: 1.0                   # Reference standard
    COSMIC_Germline: 0.9
    Commercial_Panels: 0.8
    HPO_Neoplasm: 0.7              # Lowest priority
```

### Decision Thresholds

```yaml
scoring:
  thresholds:
    score_threshold: 1.5        # Minimum score for inclusion
    watch_list_threshold: 1.0   # Threshold for emerging evidence
    min_sources: 1              # Minimum supporting sources
    max_evidence_score: 5.0     # Score normalization cap
```

### Veto System (Override)

Critical sources can override scoring thresholds:

```yaml
data_sources:
  ACMG_Incidental_Findings:
    veto:
      enabled: true
      reason: "ACMG recommended for reporting of incidental findings"
      
  Manual_Curation:
    veto:
      enabled: true
      reason: "Manually curated and reviewed by clinical experts"
```

## Data Source Configuration

### Enabling/Disabling Sources

```yaml
data_sources:
  PanelApp:
    enabled: true    # Include this source
    
  Commercial_Panels:
    enabled: false   # Skip this source
```

### Adding Custom Gene Lists

#### In-house Panels

```yaml
data_sources:
  Inhouse_Panels:
    enabled: true
    panels:
      - name: "My_Custom_Panel"
        file_path: "data/my_genes.xlsx"
        gene_column: "Gene_Symbol"
        sheet_name: "Sheet1"  # For Excel files
        evidence_score: 1.0
```

#### Manual Curation Lists

```yaml
data_sources:
  Manual_Curation:
    enabled: true
    lists:
      - name: "Literature_Review_2024"
        file_path: "data/manual/literature_genes_2024.txt"
        evidence_score: 1.2
      - name: "Expert_Panel_Recommendations"
        file_path: "data/manual/expert_panel_genes.csv"
        gene_column: "Gene"
        evidence_score: 1.1
```

### Supported File Formats

- **Excel** (.xlsx, .xls) - Specify `sheet_name` if needed
- **CSV** (.csv) - Comma-separated values
- **TSV** (.tsv) - Tab-separated values  
- **Plain text** (.txt) - One gene per line

## Output Configuration

### File Formats

```yaml
output:
  formats:
    - "excel"    # Comprehensive Excel file
    - "csv"      # CSV for programmatic use
    - "parquet"  # Efficient binary format
    - "bed"      # BED files for genome browsers
```

### BED File Generation

```yaml
output:
  bed_files:
    germline: true
    padding: 0
    exons:
      enabled: true
      canonical_transcript: true
      mane_select_transcript: true
      exon_padding: 10
```

### HTML Report

```yaml
output:
  html_report:
    enabled: true
    include_summary: true
    include_top_genes: true
    include_datatable: true
```

### Intermediate Files

For debugging and analysis:

```yaml
output:
  intermediate_files:
    enabled: true
    format: "excel"
    include_raw_data: true
    include_standardized_data: true
    include_merged_data: true
    include_scored_data: true
```

## Quality Control

```yaml
quality_control:
  require_hgnc_match: true     # Remove genes not in HGNC
  require_coordinates: true    # Remove genes without coordinates
  remove_duplicates: true     # Remove duplicates by HGNC ID
  min_gene_size: 1000         # Minimum gene size in base pairs
```

## Performance Settings

```yaml
performance:
  max_workers: 4              # Parallel workers for API calls
  batch_size: 300             # Genes per batch for lookups
  enable_caching: true        # Cache API responses
  cache_ttl: 2592000         # Cache TTL in seconds (30 days)
```

## Annotation Settings

```yaml
annotation:
  genomic_coordinates: true
  transcript_info: true
  mane_transcripts: true
  gene_descriptions: true
  transcript_padding: 25      # BP padding for transcripts
  gene_padding: 5000          # BP padding for genes
  assemblies:
    - "GRCh38"
    - "GRCh37"               # Optional legacy support
```

## Examples

### Strict Clinical Panel

For clinical use with high confidence requirements:

```yaml
scoring:
  thresholds:
    score_threshold: 2.0
    min_sources: 2

data_sources:
  Commercial_Panels:
    enabled: false            # Disable commercial panels
  
  HPO_Neoplasm:
    enabled: false            # Disable automated associations
```

### Research Panel

For research with broader gene inclusion:

```yaml
scoring:
  thresholds:
    score_threshold: 1.0
    min_sources: 1

data_sources:
  Commercial_Panels:
    enabled: true
    evidence_score: 1.0       # Increase commercial panel weight
```

### Custom Scoring Weights

Adjust for specific use cases:

```yaml
scoring:
  source_group_weights:
    Inhouse_Panels: 2.0       # Prioritize in-house expertise
    ACMG_Incidental_Findings: 1.0  # Reduce ACMG weight
```

## Validation

Always validate your configuration:

```bash
custom-panel config-check -c my_config.yml
```

## Next Steps

- **[COSMIC Setup](./cosmic_setup.md)** - Configure COSMIC access
- **[OMIM Setup](./omim_setup.md)** - Configure OMIM/HPO access
- **[Running the Pipeline](./running_pipeline.md)** - Execute with your configuration