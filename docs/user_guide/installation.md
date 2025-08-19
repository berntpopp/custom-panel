# Installation Guide

This guide will help you install Custom Panel and its dependencies using the modern uv package manager.

## Prerequisites

- **Python 3.10+** - Custom Panel requires Python 3.10 or later
- **uv** - Ultra-fast Python package manager and project manager
- **Git** - For cloning the repository
- **Make** - For running development commands (usually pre-installed)

## Quick Setup

### 1. Install uv

Install uv using the official installer:

```bash
# On macOS and Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# On Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

# Alternative: using pip
pip install uv
```

### 2. Clone and Setup

```bash
# Clone the repository
git clone https://github.com/berntpopp/custom-panel.git
cd custom-panel

# Set up development environment (recommended)
make dev-setup
```

That's it! The `make dev-setup` command will:
- Install the correct Python version with uv
- Install all dependencies (main + development + scrapers)
- Set up the project environment

## Installation Options

### Development Environment (Recommended)

```bash
make dev-setup    # Install everything for development
```

### Production Environment

```bash
make install      # Install only main dependencies
```

### With Scrapers (for web scraping features)

```bash
make install-scrapers    # Install with web scraping dependencies
```

## Alternative Installation Methods

### Using uv directly

```bash
# Install all dependencies
uv sync

# Install with scrapers group
uv sync --group scrapers

# Install development dependencies
uv sync --group dev
```

### Using pip (not recommended)

```bash
# Install in development mode (after cloning)
pip install -e .
```

## Verify Installation

### Check Installation

```bash
# Using Make (recommended)
make run-help

# Or directly
uv run custom-panel --help
```

You should see the help message with available commands.

### Configuration Check

```bash
# Using Make
make config-check

# Or directly
uv run custom-panel config-check
```

## Development Workflow

### Available Make Commands

View all available commands:

```bash
make help
```

Common development tasks:

```bash
make quality      # Run linting and type checking
make test         # Run tests  
make test-cov     # Run tests with coverage
make docs-serve   # Serve documentation locally
make clean        # Clean build artifacts
```

### Running the Application

```bash
# Using Make (recommended)
make run ARGS="run --output-dir results"

# Or directly with uv
uv run custom-panel run --output-dir results
```

## Next Steps

- **[Running the Pipeline](./running_pipeline.md)** - Learn basic usage
- **[Configuration](./configuration.md)** - Customize data sources and scoring  
- **[COSMIC Setup](./cosmic_setup.md)** - Configure COSMIC access (optional)
- **[OMIM Setup](./omim_setup.md)** - Configure OMIM access (optional)

## Troubleshooting

### Common Issues

**Command not found**: Ensure uv is in your PATH after installation.

```bash
# Check uv installation
uv --version
```

**Python version errors**: uv will automatically install the correct Python version.

**Permission errors**: Make sure you have write permissions to the project directory.

**Import errors**: Run `make install-dev` to ensure all dependencies are installed.

### Getting Help

If you encounter issues:

1. Check the [Troubleshooting Guide](../troubleshooting.md)
2. Review [GitHub Issues](https://github.com/berntpopp/custom-panel/issues)  
3. Create a new issue with detailed error messages

## Advantages of uv

- **10-100x faster** dependency resolution compared to pip/Poetry
- **Automatic Python management** - installs correct Python version
- **Universal lock files** - cross-platform compatibility
- **No virtual environment management** - handled automatically
- **Better caching** - faster subsequent installs
- **Standard Python packaging** - follows PEP standards