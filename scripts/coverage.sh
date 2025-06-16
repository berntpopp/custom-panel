#!/bin/bash
# Coverage reporting script for custom-panel
# Usage: ./scripts/coverage.sh [html|xml|term]

set -e

# Default to terminal output
OUTPUT_FORMAT="${1:-term}"

echo "🧪 Running tests with coverage..."

case "$OUTPUT_FORMAT" in
    "html")
        echo "📊 Generating HTML coverage report..."
        poetry run pytest --cov=custom_panel --cov-report=html
        echo "✅ HTML coverage report generated in htmlcov/"
        echo "📖 Open htmlcov/index.html in your browser to view the report"
        ;;
    "xml")
        echo "📊 Generating XML coverage report..."
        poetry run pytest --cov=custom_panel --cov-report=xml
        echo "✅ XML coverage report generated as coverage.xml"
        ;;
    "term")
        echo "📊 Generating terminal coverage report..."
        poetry run pytest --cov=custom_panel --cov-report=term-missing
        ;;
    *)
        echo "❌ Invalid format: $OUTPUT_FORMAT"
        echo "Usage: $0 [html|xml|term]"
        exit 1
        ;;
esac

echo "🎉 Coverage analysis complete!"