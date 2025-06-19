# Data Sources Overview

Custom Panel aggregates gene information from multiple trusted sources, each selected for specific strengths in clinical genomics. This page explains each data source, its purpose, and the rationale for inclusion.

## Source Categories

Data sources are organized into three main categories:

1. **Clinical Guidelines** - Evidence-based recommendations from professional organizations
2. **Community Curated** - Expert-reviewed databases with broad community input
3. **Automated/Commercial** - Large-scale databases and commercial test offerings

## Clinical Guidelines Sources

### ACMG Incidental Findings (SF v3.2)

**Purpose**: Genes recommended for reporting as incidental findings in clinical sequencing

**Rationale**: 
- Evidence-based recommendations from the American College of Medical Genetics
- Represents genes with actionable clinical interventions
- Updated regularly based on clinical evidence
- High clinical utility and established management guidelines

**Data Type**: 
- Fixed list of ~95 genes
- Categories: Hereditary cancer, cardiovascular, metabolic, etc.

**Veto Power**: Yes - These genes are always included regardless of other scoring

### ClinGen Gene-Disease Validity

**Purpose**: Curated gene-disease relationships with clinical validity classifications

**Rationale**:
- Rigorous curation process by clinical experts
- Standardized classification system (Definitive/Strong/Moderate/Limited)
- Focus on clinical validity for genetic testing
- Transparent evidence evaluation

**Data Type**:
- Gene-disease pairs with validity classifications
- Filtered for cancer/neoplasm phenotypes
- Updated via live API

### Manual Curation Lists

**Purpose**: Expert-curated gene lists from clinical laboratories or research groups

**Rationale**:
- Incorporates local expertise and clinical experience
- Allows inclusion of emerging evidence
- Flexible for specific clinical contexts
- Can include genes not yet in public databases

**Data Type**:
- Custom Excel/CSV/text files
- User-defined evidence scores

**Veto Power**: Yes - Manually curated genes override scoring thresholds

## Community Curated Sources

### PanelApp (Genomics England)

**Purpose**: Community-curated gene panels with evidence-based ratings

**Rationale**:
- Transparent review process by multiple experts
- Traffic light system (Green/Amber/Red) indicates confidence
- Widely used in clinical practice
- Regular updates based on new evidence

**Data Type**:
- Multiple cancer-related panels
- Evidence levels: Green (3), Amber (2), Red (1)
- Includes MOI and phenotype information

**Configuration**:
```yaml
data_sources:
  PanelApp:
    panels:
      - id: 1    # Childhood solid tumours
      - id: 2    # Adult solid tumours  
      - id: 3    # Haematological malignancies
```

### In-house Panels

**Purpose**: Gene panels currently used in clinical practice by the laboratory

**Rationale**:
- Reflects current clinical practice
- Validated for local patient populations
- Incorporates institutional expertise
- Ensures continuity of care

**Data Type**:
- Excel/CSV files with gene lists
- Can include additional metadata

### TheGenCC Database

**Purpose**: Gene-disease associations from multiple curation efforts

**Rationale**:
- Aggregates data from multiple expert curation groups
- Standardized evidence classifications
- Broader coverage than single databases
- Quality control through member organizations

**Data Type**:
- Gene-disease associations with classifications
- Filtered for cancer-related diseases

## Scientific Databases

### COSMIC Cancer Gene Census

**Purpose**: Catalog of genes with mutations causally implicated in cancer

**Rationale**:
- Comprehensive cancer gene resource
- Tier system indicates evidence strength
- Includes both somatic and germline mutations
- Widely cited in cancer research

**Data Type**:
- Genes with cancer driver mutations
- Tier 1: Well-documented cancer genes
- Tier 2: Emerging evidence

**Setup Required**: [COSMIC credentials](./cosmic_setup.md)

### HPO/OMIM Neoplasm Genes

**Purpose**: Genes associated with neoplasm phenotypes in Human Phenotype Ontology

**Rationale**:
- Comprehensive phenotype-gene relationships
- Hierarchical ontology captures related concepts
- Links to OMIM for detailed clinical information
- Useful for discovering less common associations

**Data Type**:
- Genes linked to HP:0002664 (Neoplasm) and descendants
- Requires OMIM genemap2.txt file

**Setup Required**: [OMIM access token](./omim_setup.md)

## Commercial Sources

### Commercial Gene Panels

**Purpose**: Genes included in commercial hereditary cancer test panels

**Rationale**:
- Reflects market adoption and clinical demand
- Indicates genes considered actionable by multiple labs
- Helps identify consensus across providers
- Useful for comparison with academic sources

**Included Providers** (14 total):
- Myriad Genetics - myRisk
- Blueprint Genetics - Multiple cancer panels
- Invitae - Comprehensive cancer panel
- GeneDx - Hereditary cancer panel
- Fulgent Genetics - Comprehensive panel
- And 9 others...

**Data Collection**: [Web scraping system](./web_scraping.md)

## Source Selection Criteria

Sources were selected based on:

1. **Clinical Validity** - Evidence for gene-disease relationships
2. **Update Frequency** - Regular updates with new evidence
3. **Transparency** - Clear methodology and evidence criteria
4. **Accessibility** - Available data in usable formats
5. **Complementarity** - Different perspectives and evidence types

## Source Reliability Hierarchy

The scoring system reflects source reliability:

```
Highest Trust (1.5x weight)
├── ACMG Incidental Findings

High Trust (1.2x weight)  
├── ClinGen
├── TheGenCC
├── In-house Panels
└── Manual Curation

Standard Trust (1.0x weight)
└── PanelApp

Moderate Trust (0.8-0.9x weight)
├── COSMIC (0.9x)
└── Commercial Panels (0.8x)

Lower Trust (0.7x weight)
└── HPO/OMIM Automated
```

## Veto Functionality

The veto system allows critical sources to override the normal scoring threshold:

### How It Works

1. **Normal Scoring**: Genes must meet the score threshold (default 1.5)
2. **Veto Override**: Genes from veto sources are included regardless of score
3. **Clinical Priority**: Ensures critical genes aren't excluded

### Sources with Veto Power

- **ACMG Incidental Findings**: Clinical guidelines override scoring
- **Manual Curation**: Expert judgment overrides algorithm

### Configuration

```yaml
data_sources:
  ACMG_Incidental_Findings:
    veto:
      enabled: true
      reason: "ACMG recommended for reporting"
      
  Manual_Curation:
    veto:
      enabled: true  
      reason: "Expert clinical review"
```

### Use Cases

- **Emerging Genes**: Include genes with limited evidence but clinical importance
- **Local Practice**: Ensure continuity with current clinical offerings
- **Guidelines Compliance**: Meet professional recommendations

## Adding New Sources

To add a new data source:

1. **Evaluate** against selection criteria
2. **Implement** data fetcher in `custom_panel/sources/`
3. **Configure** in `default_config.yml`
4. **Document** rationale and data type
5. **Test** integration with existing sources

## Next Steps

- [Web Scraping Guide](./web_scraping.md) - Commercial panel data collection
- [Scoring System](./scoring_system.md) - How sources contribute to final scores
- [Configuration](./configuration.md) - Customize source settings