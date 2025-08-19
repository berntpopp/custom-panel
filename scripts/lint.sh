#!/bin/bash
# Development linting script for custom-panel

echo "🧹 Running code quality checks..."

echo "📦 Checking import sorting..."
uv run ruff check . --select I001 --fix

echo "🔍 Running Ruff linter..."
uv run ruff check .

echo "🎨 Running Ruff formatter..."
uv run ruff format .

echo "🔬 Running MyPy type checker..."
uv run mypy custom_panel/

echo "✅ Code quality checks complete!"