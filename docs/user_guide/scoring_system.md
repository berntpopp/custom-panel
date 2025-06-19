# Gene Scoring System

Custom Panel implements a sophisticated scoring algorithm that balances clinical evidence, source reliability, and expert judgment to determine which genes to include in the final panel.

## Overview

The scoring system evaluates each gene across multiple dimensions:

1. **Evidence Quality** - How strong is the gene-disease association?
2. **Source Reliability** - How trustworthy is the data source?
3. **Source Consensus** - How many independent sources support the gene?
4. **Clinical Priority** - Are there clinical guidelines or expert overrides?

## Core Scoring Formula

```
Final Score = Σ(source_evidence_score × internal_confidence × source_group_weight)
```

Where:
- **source_evidence_score**: Base score from the data source (0.0-1.5)
- **internal_confidence**: Confidence based on supporting evidence (0.0-1.0)
- **source_group_weight**: Final multiplier for source priority (0.7-1.5)

## Component Breakdown

### 1. Source Evidence Scores

Each source has a base evidence score reflecting its clinical reliability:

| Source | Base Score | Rationale |
|--------|------------|-----------|
| ACMG Incidental Findings | 1.5 | Clinical guidelines with actionable interventions |
| ClinGen/TheGenCC | 1.2 | Expert curation with validity classifications |
| In-house Panels | 1.2 | Validated in clinical practice |
| Manual Curation | 1.2 | Expert reviewed for specific context |
| PanelApp | 1.0 | Community consensus standard |
| COSMIC Germline | 0.9 | Established cancer gene catalog |
| Commercial Panels | 0.8 | Market consensus, variable validation |
| HPO Neoplasm | 0.7 | Automated associations, broader coverage |

### 2. Classification Multipliers

Sources with classification systems apply additional multipliers:

#### ClinGen Gene-Disease Validity
```yaml
"Definitive": ×1.5  # Replicated evidence, expert consensus
"Strong": ×1.2      # Compelling evidence, multiple studies
"Moderate": ×1.0    # Several studies, moderate evidence
"Limited": ×0.3     # Few studies, emerging evidence
"Disputed": ×0.3    # Conflicting evidence
"Refuted": ×0.0     # Excluded from scoring
```

#### PanelApp Evidence Levels
```yaml
Green (Level 3): ×1.0   # High confidence, clinical use
Amber (Level 2): ×0.5   # Moderate confidence, evaluation
Red (Level 1): ×0.0     # Low confidence, not recommended
```

#### COSMIC Tier System
```yaml
Tier 1: ×1.0  # Well-documented cancer genes
Tier 2: ×0.8  # Strong evidence, fewer studies
Other: ×0.4   # Limited classification
```

### 3. Internal Confidence Scoring

Confidence increases with supporting evidence using normalization functions:

#### Logistic Normalization
Used for sources where confidence plateaus:

```
confidence = 1 / (1 + e^(-k×(count-x0)))
```

Where:
- **k**: Steepness of curve (how quickly confidence increases)
- **x0**: Midpoint (count at 50% confidence)

**Examples**:
- Commercial Panels: k=0.25, x0=5 (50% confidence at 5 panels)
- PanelApp: k=0.35, x0=5 (steeper curve, faster confidence)

#### Linear Normalization
Used for sources with proportional confidence:

```
confidence = min(count / max_count, 1.0)
```

**Example**:
- In-house Panels: max_count=3 (100% confidence at 3 panels)

### 4. Source Group Weights

Final multipliers reflect clinical decision-making priorities:

```yaml
ACMG_Incidental_Findings: 1.5  # Highest - clinical guidelines
ClinGen: 1.2                    # Expert validity assessment
TheGenCC: 1.2                   # Multi-group curation
Inhouse_Panels: 1.2            # Local clinical validation
Manual_Curation: 1.1           # Expert judgment
PanelApp: 1.0                  # Reference standard
COSMIC_Germline: 0.9           # Research database
Commercial_Panels: 0.8         # Market consensus
HPO_Neoplasm: 0.7             # Automated associations
```

## Scoring Examples

### Example 1: Well-Supported Gene (BRCA1)

```
Sources:
- ACMG: 1.5 × 1.0 × 1.5 = 2.25
- PanelApp (Green): 1.0 × 1.0 × 1.0 = 1.0  
- ClinGen (Definitive): 1.2 × 1.5 × 1.2 = 2.16
- Commercial (8 panels): 0.8 × 0.67 × 0.8 = 0.43
- In-house: 1.2 × 1.0 × 1.2 = 1.44

Total Score: 7.28 → INCLUDED (well above 1.5 threshold)
```

### Example 2: Emerging Gene

```
Sources:
- Manual Curation: 1.2 × 1.0 × 1.1 = 1.32
- Commercial (2 panels): 0.8 × 0.12 × 0.8 = 0.08
- HPO: 0.7 × 1.0 × 0.7 = 0.49

Total Score: 1.89 → INCLUDED (above 1.5 threshold)
```

### Example 3: Low Evidence Gene

```
Sources:
- Commercial (1 panel): 0.8 × 0.02 × 0.8 = 0.01
- HPO: 0.7 × 1.0 × 0.7 = 0.49

Total Score: 0.50 → EXCLUDED (below 1.5 threshold)
```

## Decision Thresholds

### Primary Threshold

```yaml
score_threshold: 1.5  # Minimum score for inclusion
```

Genes must achieve this score to be included in the panel.

### Watch List Threshold

```yaml
watch_list_threshold: 1.0  # Emerging evidence tracking
```

Genes between 1.0-1.5 can be tracked for future consideration.

### Additional Criteria

```yaml
min_sources: 1        # At least one source required
max_evidence_score: 5.0  # Cap to prevent score inflation
```

## Veto System

The veto system allows critical sources to override scoring:

### How Veto Works

```python
if gene in veto_sources:
    include_gene = True  # Bypass score check
else:
    include_gene = (score >= threshold)
```

### Sources with Veto Power

**ACMG Incidental Findings**
- **Reason**: Professional guidelines for reporting
- **Use Case**: Ensure guideline compliance
- **Example**: *TP53* included even if only in ACMG

**Manual Curation**
- **Reason**: Expert clinical judgment
- **Use Case**: Include emerging or locally important genes
- **Example**: Novel gene from recent literature

### Veto Configuration

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

### Veto Statistics

The system tracks veto usage:

```json
{
  "veto_stats": {
    "total_vetoed": 15,
    "by_source": {
      "ACMG_Incidental_Findings": 12,
      "Manual_Curation": 3
    }
  }
}
```

## Score Calculation Process

### Step 1: Source Aggregation

For each gene, collect evidence from all sources:

```python
gene_sources = {
    "BRCA1": [
        {"source": "ACMG", "evidence_score": 1.5},
        {"source": "PanelApp", "evidence_score": 1.0, "classification": "Green"},
        {"source": "Commercial", "panel_count": 8}
    ]
}
```

### Step 2: Apply Classifications

Multiply base scores by classification factors:

```python
if source == "ClinGen":
    score = base_score * classification_multipliers[classification]
```

### Step 3: Calculate Confidence

Apply normalization based on source type:

```python
if normalization == "logistic":
    confidence = 1 / (1 + exp(-k * (count - x0)))
elif normalization == "linear":
    confidence = min(count / max_count, 1.0)
```

### Step 4: Apply Group Weights

Multiply by final source group weight:

```python
final_score = evidence_score * confidence * group_weight
```

### Step 5: Sum and Decide

Sum all source contributions and check thresholds:

```python
total_score = sum(source_scores)
include = (total_score >= threshold) or (gene in veto_genes)
```

## Customizing the Scoring

### Adjusting Thresholds

For stricter panels:
```yaml
scoring:
  thresholds:
    score_threshold: 2.0    # Higher bar
    min_sources: 2          # Require multiple sources
```

For broader inclusion:
```yaml
scoring:
  thresholds:
    score_threshold: 1.0    # Lower threshold
    watch_list_threshold: 0.5
```

### Modifying Weights

Prioritize specific sources:
```yaml
scoring:
  source_group_weights:
    Inhouse_Panels: 2.0     # Double weight for local panels
    Commercial_Panels: 0.5  # Reduce commercial influence
```

### Changing Classifications

Adjust classification multipliers:
```yaml
data_sources:
  ClinGen:
    classification_scores:
      "Definitive": 2.0     # Increase definitive weight
      "Limited": 0.1        # Decrease limited evidence
```

## Scoring Outputs

### Summary Statistics

```json
{
  "scoring_summary": {
    "total_unique_genes": 450,
    "genes_above_threshold": 285,
    "genes_vetoed": 15,
    "final_panel_size": 300,
    "score_distribution": {
      "0-1": 45,
      "1-2": 120,
      "2-3": 180,
      "3+": 105
    }
  }
}
```

### Gene-Level Details

Each gene includes:
- `score`: Final calculated score
- `source_count`: Number of supporting sources
- `source_details`: Breakdown by source
- `veto_status`: Whether veto was applied
- `include_decision`: Final inclusion status

## Quality Assurance

### Score Validation

- **Consistency Checks**: Ensure scores match source data
- **Range Validation**: Verify scores within expected bounds
- **Audit Trail**: Track all scoring decisions

### Regular Review

- **Threshold Tuning**: Adjust based on panel performance
- **Weight Optimization**: Refine based on clinical feedback
- **Source Evaluation**: Monitor source quality over time

## Next Steps

- [Configuration Guide](./configuration.md) - Customize scoring parameters
- [Data Sources](./data_sources.md) - Understand source contributions
- [Running Pipeline](./running_pipeline.md) - Execute with custom scoring