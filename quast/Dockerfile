FROM python:3.11-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/quast" \
    org.opencontainers.image.description="python3.11-slim container with quast" \
    org.opencontainers.image.version="5.3.0"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*


# Set working directory
WORKDIR /opt

# Copy requirements first for better layer caching
RUN pip install --no-cache-dir quast

# Verify installation
RUN command -v quast &>/dev/null || exit 1

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]