FROM python:3.11-slim AS base

LABEL maintainer="Ugur Tuna"
LABEL description="APOE Genotyping Toolkit — call, feasibility, stratify"

WORKDIR /app

COPY pyproject.toml README.md ./
COPY src/ src/

RUN pip install --no-cache-dir .

ENTRYPOINT ["apoe-toolkit"]
CMD ["--help"]
