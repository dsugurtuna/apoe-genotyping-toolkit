.PHONY: help install dev test lint typecheck fmt clean docker-build docker-run

PYTHON ?= python3

help:  ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package
	$(PYTHON) -m pip install .

dev:  ## Install with development dependencies
	$(PYTHON) -m pip install -e ".[dev]"

test:  ## Run the test suite
	$(PYTHON) -m pytest tests/ -v --tb=short

lint:  ## Run ruff linter
	$(PYTHON) -m ruff check src/ tests/

typecheck:  ## Run mypy type checker
	$(PYTHON) -m mypy src/

fmt:  ## Auto-format code with ruff
	$(PYTHON) -m ruff format src/ tests/
	$(PYTHON) -m ruff check --fix src/ tests/

clean:  ## Remove build artefacts
	rm -rf build/ dist/ *.egg-info src/*.egg-info .mypy_cache .pytest_cache .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} +

docker-build:  ## Build Docker image
	docker build -t apoe-toolkit:latest .

docker-run:  ## Run the toolkit in Docker (pass ARGS="call --help")
	docker run --rm apoe-toolkit:latest $(ARGS)
