# COSMIC Authentication Setup

This guide explains how to configure COSMIC Cancer Gene Census authentication for the custom-panel tool.

## Overview

COSMIC (Catalogue of Somatic Mutations in Cancer) requires authentication to access their Cancer Gene Census data. The custom-panel tool supports automatic login and authenticated downloads for germline cancer gene data.

## Prerequisites

1. **COSMIC Account**: You need a free COSMIC account from https://cancer.sanger.ac.uk/cosmic/register
2. **Valid Credentials**: Email and password for your COSMIC account

## Configuration

### Step 1: Create Local Configuration File

Create a `config.local.yml` file in your project root (this file should NOT be committed to version control):

```yaml
data_sources:
  cosmic_germline:
    enabled: true
    email: your-cosmic-email@example.com
    password: your-cosmic-password
    cache_dir: .cache/cosmic
    cache_expiry_days: 30
    germline_scoring:
      enabled: true
      tier_weights:
        "Tier 1": 1.0
        "Tier 2": 0.8
        "": 0.4
```

### Step 2: Replace Credentials

Replace the placeholder values:
- `your-cosmic-email@example.com` → Your actual COSMIC account email
- `your-cosmic-password` → Your actual COSMIC account password

### Step 3: Secure the Configuration File

Add `config.local.yml` to your `.gitignore` file to prevent accidentally committing credentials:

```bash
echo "config.local.yml" >> .gitignore
```

## Configuration Options

### Required Settings

- **`email`**: Your COSMIC account email address
- **`password`**: Your COSMIC account password

### Optional Settings

- **`enabled`**: Enable/disable COSMIC data source (default: `false`)
- **`cache_dir`**: Directory to store downloaded files (default: `.cache/cosmic`)
- **`cache_expiry_days`**: Days before cache expires (default: `30`)

### Scoring Configuration

#### Germline Scoring
- **`germline_scoring.enabled`**: Enable germline gene scoring
- **`germline_scoring.tier_weights`**: Evidence weights by COSMIC tier


### Tier Weights

COSMIC genes are classified into tiers based on evidence strength:
- **Tier 1**: Highest confidence cancer genes (weight: 1.0)
- **Tier 2**: Moderate confidence cancer genes (weight: 0.8)  
- **""** (empty): Unknown or unclassified genes (weight: 0.4)

## Usage

### Fetch COSMIC Data Only

```bash
# Option 1: Automatic config loading (recommended)
custom-panel fetch cosmic

# Option 2: Explicitly specify local config file
custom-panel fetch cosmic -c config.local.yml

# Or if using Python directly:
python -c "from custom_panel.cli import app; app()" fetch cosmic
```

### Run Full Pipeline (includes COSMIC if enabled)

```bash
# Option 1: Automatic config loading (recommended)
custom-panel run

# Option 2: Explicitly specify local config file  
custom-panel run -c config.local.yml

# Or if using Python directly:
python -c "from custom_panel.cli import app; app()" run
```

**Note**: The tool always loads the default configuration first, then applies overrides from either the automatically detected `config.local.yml` or the file specified with `-c`.

### Validate Configuration

```bash
python -c "from custom_panel.cli import app; app()" config-check
```

## Authentication Flow

1. **Login**: Tool authenticates with COSMIC using your credentials
2. **Session**: Maintains authenticated session with proper cookies/headers
3. **Download URL**: Requests authenticated download URL for Cancer Gene Census
4. **Download**: Downloads file using authenticated URL
5. **Cache**: Stores file locally for future use (respects cache expiry)

## Troubleshooting

### Authentication Errors

If you see authentication errors:

1. **Verify Credentials**: Double-check your email and password
2. **Check Account**: Ensure your COSMIC account is active
3. **Network Issues**: Check internet connectivity
4. **COSMIC Status**: Verify COSMIC website is accessible

### Configuration Validation

Common configuration errors:

```bash
# Missing credentials
COSMIC requires either credentials (email/password) or legacy census_url

# Invalid email format  
COSMIC email appears to be invalid

# Missing password
COSMIC email provided but password missing
```

### Debug Logging

Enable debug logging for detailed authentication information:

```bash
python -c "from custom_panel.cli import app; app()" run --log-level DEBUG
```

## Data Processing

The tool processes COSMIC data focusing on germline genes associated with inherited cancer predisposition. Each gene receives an evidence score based on its COSMIC tier classification and your configured weights.

## Security Notes

- **Never commit** `config.local.yml` to version control
- **Restrict file permissions** on the config file: `chmod 600 config.local.yml`
- **Regularly rotate** your COSMIC password
- **Use strong passwords** for your COSMIC account

## Legacy Support

For backward compatibility, you can still use direct URLs:

```yaml
data_sources:
  cosmic_germline:
    enabled: true
    census_url: https://direct-url-to-cosmic-file.csv
```

However, authenticated access is recommended for reliability and compliance with COSMIC's terms of service.

## Example Complete Configuration

```yaml
# config.local.yml
data_sources:
  cosmic_germline:
    enabled: true
    email: researcher@university.edu
    password: MySecurePassword123!
    cache_dir: .cache/cosmic
    cache_expiry_days: 30
    germline_scoring:
      enabled: true
      tier_weights:
        "Tier 1": 1.0
        "Tier 2": 0.8
        "": 0.4

  # Other data sources...
  panelapp:
    enabled: true
    # ... other configs
```