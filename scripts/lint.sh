#!/bin/bash
# Development linting script for custom-panel

echo "🧹 Running code quality checks..."

echo "📦 Checking import sorting..."
poetry run ruff check . --select I001 --fix

echo "🔍 Running Ruff linter..."
poetry run ruff check .

echo "🎨 Running Ruff formatter..."
poetry run ruff format .

echo "🔬 Running MyPy type checker..."
poetry run mypy custom_panel/

echo "✅ Code quality checks complete!"