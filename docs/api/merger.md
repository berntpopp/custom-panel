# Panel Merger Engine

The merger engine combines data from multiple sources and applies the scoring algorithm to determine final gene inclusion.

## Overview

The `PanelMerger` class implements the core scoring logic:

- **Data standardization** - Ensures consistent data formats across sources
- **Evidence aggregation** - Combines evidence from multiple sources per gene
- **Scoring algorithm** - Applies configurable weights and thresholds
- **Decision logic** - Determines final gene inclusion with veto capabilities
- **Quality metrics** - Provides detailed scoring statistics

## Scoring System

The merger implements a sophisticated scoring system that balances:

1. **Source reliability** - Evidence scores based on clinical trust
2. **Classification quality** - Multipliers for evidence classifications  
3. **Source diversity** - Confidence scores based on supporting sources
4. **Clinical priorities** - Final weights and veto capabilities

## API Reference

::: custom_panel.engine.merger
    options:
      show_root_heading: true
      show_source: false
      heading_level: 2