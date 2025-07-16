# Troubleshooting

This guide helps you resolve common issues when using the custom-panel tool.

## Common Issues

### Installation Problems

#### Poetry Installation Issues

**Problem**: Poetry installation fails or commands don't work

**Solution**:
```bash
# Install Poetry using the official installer
curl -sSL https://install.python-poetry.org | python3 -

# Or using pip
pip install poetry

# Verify installation
poetry --version
```

#### Python Version Issues

**Problem**: Python version compatibility errors

**Solution**:
- Ensure you're using Python 3.10 or later
- Use pyenv or conda to manage Python versions:

```bash
# Check Python version
python --version

# Using pyenv
pyenv install 3.10.0
pyenv local 3.10.0

# Using conda
conda create -n custom-panel python=3.10
conda activate custom-panel
```

#### Dependency Conflicts

**Problem**: Package dependency conflicts during installation

**Solution**:
```bash
# Clear poetry cache
poetry cache clear --all .

# Install with fresh resolver
poetry install --no-cache

# Install specific groups
poetry install --with dev,scrapers
```

### Configuration Issues

#### Missing Configuration Files

**Problem**: Configuration file not found errors

**Solution**:
- Ensure `config/default_config.yml` exists in the project
- Create a local configuration file:

```bash
cp config/default_config.yml config.local.yml
```

#### Invalid Configuration Format

**Problem**: YAML parsing errors

**Solution**:
- Check YAML syntax using an online validator
- Ensure proper indentation (spaces, not tabs)
- Validate configuration:

```bash
custom-panel config-check
```

#### API Key Issues

**Problem**: API authentication failures

**Solution**:
- Set up API keys in `config.local.yml`:

```yaml
# Example for COSMIC
data_sources:
  cosmic_germline:
    api_key: "your-cosmic-api-key"
    
# Example for OMIM
data_sources:
  hpo_neoplasm:
    omim_genemap2_url: "https://data.omim.org/downloads/YOUR_TOKEN/genemap2.txt"
```

### Runtime Issues

#### Network Connection Problems

**Problem**: API calls fail due to network issues

**Solution**:
- Check internet connectivity
- Verify firewall settings
- Increase timeout values in configuration:

```yaml
apis:
  ensembl:
    timeout: 60
    max_retries: 5
```

#### Memory Issues

**Problem**: Out of memory errors during processing

**Solution**:
- Reduce batch sizes:

```yaml
performance:
  batch_size: 100
  max_workers: 2
```

- Process smaller datasets
- Increase system memory or use a machine with more RAM

#### File Permission Issues

**Problem**: Cannot write to output directory

**Solution**:
```bash
# Check permissions
ls -la results/

# Fix permissions
chmod 755 results/
sudo chown -R $USER:$USER results/
```

### Data Source Issues

#### Source Unavailable

**Problem**: External data sources return errors

**Solution**:
- Check if the source is temporarily down
- Disable problematic sources:

```yaml
data_sources:
  problematic_source:
    enabled: false
```

- Use local fallback files if available

#### Empty Results

**Problem**: No data returned from sources

**Solution**:
- Check source configuration
- Verify API endpoints are correct
- Enable debug logging:

```bash
custom-panel run --log-level DEBUG
```

#### Data Format Changes

**Problem**: Source data format has changed

**Solution**:
- Check for tool updates
- Report the issue on GitHub
- Temporarily disable the affected source

### SNP Processing Issues

#### Harmonization Failures

**Problem**: SNP coordinate harmonization fails

**Solution**:
- Check Ensembl API availability
- Disable harmonization temporarily:

```yaml
snp_processing:
  harmonization:
    enabled: false
```

#### Large SNP Datasets

**Problem**: SNP processing is slow or fails

**Solution**:
- Reduce batch sizes:

```yaml
snp_processing:
  harmonization:
    ensembl_batch_size: 10
```

- Enable caching:

```yaml
snp_processing:
  harmonization:
    caching:
      enabled: true
```

### Output Issues

#### HTML Report Generation Fails

**Problem**: HTML reports are not generated

**Solution**:
- Check if HTML output is enabled:

```yaml
output:
  html_report:
    enabled: true
```

- Verify JavaScript dependencies
- Check browser compatibility

#### Excel Export Issues

**Problem**: Excel files are corrupted or won't open

**Solution**:
- Install openpyxl:

```bash
poetry add openpyxl
```

- Check file permissions
- Try different Excel formats

#### BED File Issues

**Problem**: BED files are malformed

**Solution**:
- Verify coordinate data is complete
- Check chromosome naming conventions
- Validate BED format:

```bash
# Check BED file format
head -n 5 results/*/06_final_output/germline_panel.bed
```

## Getting Help

### Enable Debug Logging

For detailed troubleshooting information:

```bash
custom-panel run --log-level DEBUG --output-dir debug_results
```

### Check Log Files

Log files are saved in the output directory:

```
results/
└── run_YYYYMMDD_HHMMSS/
    └── logs/
        ├── pipeline.log
        ├── sources.log
        └── annotator.log
```

### System Information

When reporting issues, include:

```bash
# System info
python --version
poetry --version
custom-panel --version  # If available

# Package info
poetry show
```

### Report Issues

If you can't resolve the issue:

1. Check existing issues on GitHub
2. Create a new issue with:
   - Error message
   - System information
   - Configuration file (redacted)
   - Log files (if relevant)

### Community Support

- GitHub Discussions for questions
- GitHub Issues for bug reports
- Documentation for usage guidance

## Performance Optimization

### Speed up Processing

```yaml
performance:
  max_workers: 6          # Increase parallel processing
  batch_size: 500         # Larger batches for APIs
  enable_caching: true    # Enable all caching
```

### Reduce Memory Usage

```yaml
performance:
  max_workers: 2          # Reduce parallel processing
  batch_size: 50          # Smaller batches
  
quality_control:
  min_gene_size: 5000     # Filter small genes
```

### Optimize Output

```yaml
output:
  formats:
    - "excel"             # Only generate needed formats
    - "bed"
  
  intermediate_files:
    enabled: false        # Disable intermediate files
```

## Common Error Messages

### "No data fetched from any source"

**Cause**: All data sources are disabled or failing

**Solution**:
- Enable at least one data source
- Check network connectivity
- Verify API keys and configurations

### "Failed to standardize gene symbols"

**Cause**: HGNC API is unavailable or rate-limited

**Solution**:
- Check internet connection
- Reduce batch sizes
- Enable caching for HGNC calls

### "No genes in master list after merging"

**Cause**: Score threshold is too high or no genes meet criteria

**Solution**:
- Lower score threshold:

```yaml
scoring:
  thresholds:
    score_threshold: 1.0
```

### "Annotation failed"

**Cause**: Ensembl API issues or invalid gene symbols

**Solution**:
- Check Ensembl API status
- Verify gene symbol standardization
- Enable annotation fallbacks

## Best Practices

### Configuration Management

1. Use `config.local.yml` for environment-specific settings
2. Keep sensitive data (API keys) in local config only
3. Version control your configuration changes
4. Test configuration with `config-check` before running

### Error Handling

1. Always check logs for detailed error information
2. Use debug mode for troubleshooting
3. Test with small datasets first
4. Keep backups of working configurations

### Performance

1. Start with default settings
2. Adjust batch sizes based on your system
3. Enable caching for repeated runs
4. Monitor memory usage with large datasets

### Updates

1. Check for tool updates regularly
2. Review changelog for breaking changes
3. Test updates with known configurations
4. Report issues encountered with new versions