.PHONY: help install install-dev install-scrapers clean lint format typecheck test test-cov docs docs-serve build run sync lock
.DEFAULT_GOAL := help

# Colors for output
GREEN := \033[32m
YELLOW := \033[33m
BLUE := \033[34m
RED := \033[31m
NC := \033[0m # No Color

help: ## Show this help message
	@echo "$(BLUE)Custom Panel - Development Commands$(NC)"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "$(GREEN)%-20s$(NC) %s\n", $$1, $$2}'

# Environment Management
install: ## Install main dependencies
	@echo "$(YELLOW)📦 Installing main dependencies...$(NC)"
	uv sync --no-dev

install-dev: ## Install all dependencies including dev
	@echo "$(YELLOW)📦 Installing all dependencies (including dev)...$(NC)"
	uv sync

install-scrapers: ## Install all dependencies (scrapers now included in main)
	@echo "$(YELLOW)📦 Installing all dependencies (scrapers included in main)...$(NC)"
	uv sync

sync: ## Sync dependencies with lock file
	@echo "$(YELLOW)🔄 Syncing dependencies...$(NC)"
	uv sync

lock: ## Update lock file
	@echo "$(YELLOW)🔒 Updating lock file...$(NC)"
	uv lock

# Code Quality
lint: ## Run linting (ruff check)
	@echo "$(YELLOW)🔍 Running Ruff linter...$(NC)"
	uv run ruff check .

format: ## Format code (ruff format)
	@echo "$(YELLOW)🎨 Running Ruff formatter...$(NC)"
	uv run ruff format .

format-check: ## Check if code is formatted correctly
	@echo "$(YELLOW)🎨 Checking code formatting...$(NC)"
	uv run ruff format --check .

typecheck: ## Run type checking (mypy)
	@echo "$(YELLOW)🔬 Running MyPy type checker...$(NC)"
	uv run mypy custom_panel/

quality: lint format-check typecheck ## Run all code quality checks

# Testing
test: ## Run tests
	@echo "$(YELLOW)🧪 Running tests...$(NC)"
	uv run pytest -v

test-cov: ## Run tests with coverage
	@echo "$(YELLOW)🧪 Running tests with coverage...$(NC)"
	uv run pytest -v --cov=custom_panel --cov-report=term-missing

test-cov-html: ## Run tests with HTML coverage report
	@echo "$(YELLOW)🧪 Running tests with HTML coverage...$(NC)"
	uv run pytest -v --cov=custom_panel --cov-report=html
	@echo "$(GREEN)📊 Coverage report generated in htmlcov/$(NC)"

test-cov-xml: ## Run tests with XML coverage report
	@echo "$(YELLOW)🧪 Running tests with XML coverage...$(NC)"
	uv run pytest -v --cov=custom_panel --cov-report=xml
	@echo "$(GREEN)📊 Coverage report generated as coverage.xml$(NC)"

# Documentation
docs: ## Build documentation
	@echo "$(YELLOW)📚 Building documentation...$(NC)"
	uv run mkdocs build

docs-serve: ## Serve documentation locally
	@echo "$(YELLOW)📚 Serving documentation at http://127.0.0.1:8000$(NC)"
	uv run mkdocs serve

# Application Commands
run: ## Run the application (requires arguments)
	@echo "$(YELLOW)🚀 Running custom-panel...$(NC)"
	@echo "$(RED)Usage: make run ARGS='<command> <options>'$(NC)"
	@echo "$(RED)Example: make run ARGS='fetch panelapp --output-dir results'$(NC)"

run-help: ## Show application help
	@echo "$(YELLOW)📖 Showing custom-panel help...$(NC)"
	uv run custom-panel --help

config-check: ## Check configuration
	@echo "$(YELLOW)⚙️  Checking configuration...$(NC)"
	uv run custom-panel config-check

# Utility Commands
clean: ## Clean build artifacts and cache
	@echo "$(YELLOW)🧹 Cleaning build artifacts...$(NC)"
	rm -rf dist/
	rm -rf build/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf htmlcov/
	rm -f coverage.xml
	rm -f .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build: ## Build the package
	@echo "$(YELLOW)🏗️  Building package...$(NC)"
	uv build

# Development Workflow
dev-setup: install-dev ## Set up development environment
	@echo "$(GREEN)✅ Development environment ready!$(NC)"
	@echo "$(BLUE)Run 'make help' to see available commands$(NC)"

ci: quality test ## Run all CI checks locally
	@echo "$(GREEN)✅ All CI checks passed!$(NC)"

# Custom commands with arguments
ifdef ARGS
run:
	@echo "$(YELLOW)🚀 Running: custom-panel $(ARGS)$(NC)"
	uv run custom-panel $(ARGS)
endif