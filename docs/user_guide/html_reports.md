# HTML Reports User Guide

The custom-panel tool generates rich, interactive HTML reports that provide comprehensive visualizations and analysis of your gene panel data. This guide explains how to enable, generate, and use these reports.

## Overview

The HTML report is an interactive web-based dashboard that includes:

- **Gene Data Tables**: Sortable, filterable tables with all gene information
- **SNP Data Tables**: Comprehensive SNP analysis across multiple categories
- **Interactive Charts**: Visual representations of data distributions and statistics
- **Export Functionality**: Download capabilities for Excel files
- **Source Analysis**: Gene data organized by source for easy exploration
- **Statistics Dashboard**: Key metrics and summary information

## Enabling HTML Reports

HTML reports are generated automatically when enabled in your configuration. 

### Configuration Options

In your `config.local.yml` or custom configuration file:

```yaml
output:
  html:
    enabled: true  # Enable HTML report generation
```

### Default Settings

HTML reports are enabled by default in the standard configuration. You can verify this by running:

```bash
custom-panel config-check
```

## Generating Reports

HTML reports are generated as part of the standard pipeline execution:

```bash
# Generate panel with HTML report
custom-panel run --output-dir results

# Generate with debug information
custom-panel run --output-dir results --log-level DEBUG
```

The HTML report will be saved as `panel_report.html` in the final output directory:

```
results/
└── run_YYYYMMDD_HHMMSS/
    └── 06_final_output/
        └── panel_report.html  # ← Interactive HTML report
```

## Report Features

### 1. Gene Data Analysis

#### Interactive Gene Tables

- **All Genes Tab**: Complete gene list with filtering and sorting
- **Panel Genes Tab**: Genes recommended for inclusion in the panel
- **Source-Specific Tabs**: Individual data source breakdowns (PanelApp, ACMG, etc.)

#### Gene Table Columns

- **Gene Symbol**: HGNC gene symbol with external links
- **Gene Name**: Full gene description
- **Chromosome**: Chromosomal location
- **Score**: Calculated evidence score
- **Include**: Panel inclusion recommendation
- **Sources**: Contributing data sources
- **Coordinates**: Genomic coordinates (GRCh38)

#### Interactive Features

- **Search**: Real-time search across all columns
- **Sorting**: Click column headers to sort data
- **Filtering**: Use column filters for precise data selection
- **Pagination**: Navigate through large datasets
- **Export**: Download filtered data as Excel files

### 2. SNP Data Analysis

The HTML report includes comprehensive SNP analysis across multiple categories:

#### SNP Categories

- **Identity SNPs**: Forensic identification markers
- **Ethnicity SNPs**: Ancestry inference markers  
- **PRS SNPs**: Polygenic risk score variants
- **ClinVar SNPs**: Clinical significance variants
- **Manual SNPs**: Curator-defined variants

#### SNP Table Features

- **Variant Information**: rsID, coordinates, alleles
- **Clinical Data**: Significance, phenotypes, populations
- **Quality Metrics**: Validation status and confidence scores
- **Source Attribution**: Origin database and references

### 3. Interactive Charts

The report includes several interactive visualizations:

#### Score Distribution Chart

- Shows distribution of gene evidence scores
- Helps identify highly-supported genes
- Interactive hover tooltips with detailed information

#### Gene Size Distribution

- Visualizes gene length distributions
- Identifies unusually large or small genes
- Useful for panel design considerations

#### Source Count Analysis

- Shows how many sources support each gene
- Helps assess evidence strength
- Identifies consensus genes across databases

#### Transcript Size Distribution

- Displays transcript length variations
- Useful for sequencing coverage planning
- Helps identify potential technical challenges

### 4. Statistics Dashboard

The report header provides key summary statistics:

- **Total Genes**: Complete gene count
- **Panel Genes**: Recommended inclusion count
- **Average Score**: Mean evidence score
- **Source Coverage**: Number of contributing databases
- **SNP Counts**: Variants by category
- **Quality Metrics**: Data completeness indicators

### 5. Export Functionality

#### Excel Export Options

- **Gene Data**: Complete gene information with all annotations
- **SNP Data**: Full variant details with clinical information
- **Filtered Data**: Export only currently filtered/visible data
- **Source-Specific**: Export individual data source contents

#### Export Process

1. Apply desired filters to the data tables
2. Click the "Export to Excel" button
3. Select export scope (all data vs. filtered data)
4. Download automatically starts

## Browser Compatibility

The HTML reports are optimized for modern web browsers:

- **Chrome**: Fully supported (recommended)
- **Firefox**: Fully supported
- **Safari**: Fully supported
- **Edge**: Fully supported
- **Internet Explorer**: Not supported

### JavaScript Requirements

The reports require JavaScript to be enabled for full functionality:

- Interactive tables and filtering
- Chart rendering and interactions
- Export functionality
- Dynamic content loading

## Performance Considerations

### Large Datasets

For datasets with thousands of genes:

- Tables use pagination for optimal performance
- Charts sample data points for smooth rendering
- Export functions handle large datasets efficiently
- Browser memory usage is optimized

### Loading Times

- Initial page load: ~1-3 seconds for typical datasets
- Interactive features: Near-instantaneous response
- Export operations: 5-10 seconds for large datasets
- Chart rendering: ~1-2 seconds

## Customization Options

### Styling

The reports use modern CSS with:

- Responsive design for mobile devices
- High contrast for accessibility
- Professional color scheme
- Print-friendly layouts

### Data Display

Configure what data appears in reports through your configuration:

```yaml
output:
  html:
    enabled: true
    include_charts: true      # Include interactive charts
    include_statistics: true  # Include summary statistics
    max_genes_display: 10000  # Limit for performance
```

## Troubleshooting

### Common Issues

#### Report Not Generated

1. Check if HTML output is enabled in configuration
2. Verify sufficient disk space in output directory
3. Check log files for error messages
4. Ensure all required dependencies are installed

#### JavaScript Errors

1. Try refreshing the page
2. Check browser console for error messages
3. Verify JavaScript is enabled in browser settings
4. Try a different browser

#### Export Not Working

1. Check browser popup/download settings
2. Verify sufficient disk space
3. Try exporting smaller data subsets
4. Check browser compatibility

#### Charts Not Displaying

1. Ensure JavaScript is enabled
2. Check browser console for errors
3. Verify internet connection (for CDN resources)
4. Try refreshing the page

### Getting Help

If you encounter issues with HTML reports:

1. Check the log files in your output directory
2. Review the [troubleshooting guide](../troubleshooting.md)
3. Consult the [API documentation](../api/cli.md) for technical details
4. Report issues on the project GitHub repository

## Best Practices

### Report Usage

1. **Review Statistics First**: Check summary statistics before diving into detailed data
2. **Use Filters Effectively**: Apply filters to focus on specific gene sets or criteria
3. **Export Strategically**: Export filtered data rather than complete datasets when possible
4. **Monitor Performance**: For large datasets, use pagination and filtering to maintain responsiveness

### Data Interpretation

1. **Score Analysis**: Higher scores indicate stronger evidence for panel inclusion
2. **Source Consensus**: Genes supported by multiple sources are generally more reliable
3. **Quality Metrics**: Review completeness and validation status indicators
4. **Clinical Relevance**: Cross-reference with clinical guidelines and literature

### Integration Workflow

1. **Generate Reports Regularly**: Include HTML reports in your standard pipeline runs
2. **Share with Team**: HTML reports are ideal for sharing results with colleagues
3. **Archive Important Reports**: Save reports for future reference and comparison
4. **Link to Documentation**: Reference this guide when sharing reports with new users

## Sample Report

A [sample HTML report](../examples/sample_panel_report.html) is available to demonstrate all features and functionality described in this guide.