# OMIM Configuration Setup

This document explains how to configure access to OMIM (Online Mendelian Inheritance in Man) data for the HPO/OMIM neoplasm gene identification feature.

## Overview

The HPO/OMIM data source identifies cancer-related genes by:
1. Fetching all descendant terms of HP:0002664 (Neoplasm) from the HPO API
2. Downloading the phenotype.hpoa file (automatically cached)
3. Filtering phenotype annotations for neoplasm-related HPO terms
4. Cross-referencing with OMIM genemap2 data to find associated genes

## OMIM Access Requirements

OMIM requires registration and personal access tokens for downloading their genemap2.txt file, which contains gene-to-disease mappings.

### Step 1: Register with OMIM

1. Go to https://omim.org/downloads
2. Register for an account (free for academic/research use)
3. Accept the terms and conditions
4. You will receive a personal access token/URL

### Step 2: Configure Your Access Token

You have two options:

#### Option A: Local Configuration File (Recommended)

1. Copy the template configuration:
   ```bash
   cp config.local.yml.template config.local.yml
   ```

2. Edit `config.local.yml` and replace `YOUR_TOKEN` with your actual token:
   ```yaml
   data_sources:
     hpo_neoplasm:
       omim_genemap2_url: "https://data.omim.org/downloads/YOUR_ACTUAL_TOKEN/genemap2.txt"
   ```

   Example with a real token:
   ```yaml
   data_sources:
     hpo_neoplasm:
       omim_genemap2_url: "https://data.omim.org/downloads/9GJLEFvqSmWaImCijeRdVA/genemap2.txt"
   ```

#### Option B: Environment Variable

Set the URL as an environment variable:
```bash
export OMIM_GENEMAP2_URL="https://data.omim.org/downloads/YOUR_TOKEN/genemap2.txt"
```

Then in your `config.local.yml`:
```yaml
data_sources:
  hpo_neoplasm:
    omim_genemap2_url: "${OMIM_GENEMAP2_URL}"
```

#### Option C: Manual Download

If you prefer to download the file manually:

1. Download genemap2.txt from your OMIM downloads page
2. Place it in `data/omim/genemap2.txt`
3. The system will use the local file automatically

## Configuration Options

All HPO/OMIM configuration options:

```yaml
data_sources:
  hpo_neoplasm:
    enabled: true
    
    # Option 1: Direct URL with access token (recommended)
    omim_genemap2_url: "https://data.omim.org/downloads/YOUR_TOKEN/genemap2.txt"
    
    # Option 2: Local file path (fallback)
    omim_genemap2_path: "data/omim/genemap2.txt"
    
    # Cache settings
    cache_dir: ".cache/hpo"
    cache_expiry_days: 30
    
    # Scoring parameters
    base_evidence_score: 0.7
```

## Security Notes

- **Never commit your access token to version control**
- The `config.local.yml` file is automatically ignored by git
- Keep your access token private and secure
- OMIM tokens may expire and need renewal

## Usage

Once configured, use the HPO data source:

```bash
# Fetch HPO/OMIM neoplasm genes
custom-panel fetch hpo --output-dir results/hpo

# Use in full pipeline
custom-panel run --output-dir results/final
```

## Troubleshooting

### "No OMIM genemap2 URL or file path configured"
- Ensure you have created `config.local.yml` with your access token
- Check that the token URL is correct

### "Failed to download OMIM genemap2"
- Verify your access token is still valid
- Check your internet connection
- Try downloading manually from the OMIM website

### "OMIM genemap2 file not found"
- If using local file option, ensure the file exists at the specified path
- Check file permissions

### Empty results
- Verify the genemap2.txt file downloaded correctly (should be ~3MB)
- Check that the HPO API is accessible

## Expected Results

When working correctly, the HPO/OMIM data source should identify:
- ~15-50 cancer-related genes (depending on genemap2.txt version)
- Well-known cancer genes like TP53, BRCA1/2, MLH1, MSH2, APC, VHL, etc.
- Evidence scores based on number of associated HPO terms and diseases

## File Locations

- Configuration: `config.local.yml` (create from template)
- Template: `config.local.yml.template`
- Cache: `.cache/hpo/` directory
- Local files: `data/omim/` directory