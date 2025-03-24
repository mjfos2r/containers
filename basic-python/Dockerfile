FROM python:3.11-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/basic-python" \
    org.opencontainers.image.description="Basic python3.11-slim container with pigz and other utils." \
    org.opencontainers.image.version="3.11-slim"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    pigz \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

ENTRYPOINT [ "/bin/bash" ]
