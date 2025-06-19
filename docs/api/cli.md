# CLI Reference

The Custom Panel command-line interface provides powerful tools for gene panel curation and management.

## Overview

The CLI is built with Typer and provides several main commands:

- **`run`** - Execute the complete gene panel pipeline
- **`fetch`** - Fetch data from individual sources  
- **`config-check`** - Validate configuration files
- **`search-panels`** - Search available panels in PanelApp

## Main Module

::: custom_panel.cli
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2