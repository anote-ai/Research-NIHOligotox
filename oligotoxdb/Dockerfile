##############################################################################
# OligoToxDB — Reproducibility Container
#
# Provides a fully reproducible environment for:
#   - Running the data ingestion pipeline
#   - Training OligoTox-XGB and OligoTox-Transformer models
#   - Launching the OligoTox-Web portal
#   - Running the full test suite
#
# Build:  docker build -t oligotoxdb:latest .
# Run pipeline:
#   docker run -v $(pwd)/data:/data oligotoxdb python scripts/generate_synthetic_data.py
# Run portal:
#   docker run -p 7860:7860 -e OLIGOTOXDB_PATH=/data/db.duckdb oligotoxdb oligotox-portal
# Run tests:
#   docker run oligotoxdb python -m pytest tests/ -v
##############################################################################

FROM python:3.11-slim-bookworm

LABEL maintainer="oligotoxdb@anote.ai"
LABEL version="0.1.0"
LABEL description="OligoToxDB: Open oligonucleotide toxicity dataset and ML toolkit"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd -m -u 1000 oligotox
WORKDIR /app

# Install Python dependencies in stages for better layer caching
# Stage 1: core scientific stack
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
        numpy>=1.26 \
        pandas>=2.1 \
        scipy>=1.11 \
        scikit-learn>=1.4 \
        xgboost>=2.0 \
        shap>=0.44 \
        duckdb>=0.10 \
        pydantic>=2.5 \
        click>=8.1 \
        rich>=13.0 \
        tqdm>=4.66

# Stage 2: ML (optional — comment out to reduce image size)
RUN pip install --no-cache-dir \
        torch>=2.2 --index-url https://download.pytorch.org/whl/cpu \
        transformers>=4.40 \
        datasets>=2.18 \
        accelerate>=0.29 \
        GPy>=1.13

# Stage 3: Portal and visualization
RUN pip install --no-cache-dir \
        gradio>=4.26 \
        plotly>=5.20 \
        matplotlib>=3.8 \
        seaborn>=0.13

# Copy source code
COPY . .
RUN pip install --no-cache-dir -e . 2>/dev/null || true

# Set ownership
RUN chown -R oligotox:oligotox /app
USER oligotox

# Environment
ENV PYTHONPATH=/app
ENV OLIGOTOXDB_PATH=/data/oligotoxdb.duckdb
ENV OLIGOTOX_MODEL_DIR=/models/xgb

# Data and model volumes
VOLUME ["/data", "/models"]

# Health check for portal
HEALTHCHECK --interval=30s --timeout=10s --start-period=15s --retries=3 \
    CMD curl -f http://localhost:7860/ || exit 1

# Default: run tests
CMD ["python", "-m", "pytest", "tests/", "-v", "--tb=short"]

# To run the portal instead:
# docker run -p 7860:7860 oligotoxdb python portal/app.py
