# Installation Guide

This guide will help you install Custom Panel and its dependencies.

## Prerequisites

- **Python 3.10+** - Custom Panel requires Python 3.10 or later
- **Poetry** - Recommended package manager for Python projects
- **Git** - For cloning the repository

## Installation Methods

### Using Poetry (Recommended)

First, clone the repository:

```bash
git clone https://github.com/berntpopp/custom-panel.git
cd custom-panel
```

Install with Poetry:

```bash
poetry install
poetry shell
```

### Using pip

```bash
# Install in development mode (after cloning)
pip install -e .
```

### From PyPI (when published)

```bash
pip install custom-panel
```

## Verify Installation

Check that the installation was successful:

```bash
custom-panel --help
```

You should see the help message with available commands.

## Configuration Check

Verify your configuration is valid:

```bash
custom-panel config-check
```

## Next Steps

- **[Running the Pipeline](./running_pipeline.md)** - Learn basic usage
- **[Configuration](./configuration.md)** - Customize data sources and scoring
- **[COSMIC Setup](./cosmic_setup.md)** - Configure COSMIC access (optional)
- **[OMIM Setup](./omim_setup.md)** - Configure OMIM access (optional)

## Troubleshooting

### Common Issues

**ImportError**: Ensure you're using Python 3.10+ and all dependencies are installed.

**Permission errors**: Make sure you have write permissions to the output directory.

**API timeouts**: Check your internet connection and API access credentials.

### Getting Help

If you encounter issues:

1. Check the [GitHub Issues](https://github.com/berntpopp/custom-panel/issues)
2. Review the configuration documentation
3. Create a new issue with detailed error messages