FROM python:3.12-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers" \
    org.opencontainers.image.description="Python3.12-slim container with Bam2Plot and other plotting tools included" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="your.email@example.com"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    zlib1g-dev \ 
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# install mosdepth 0.3.11
RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.11/mosdepth \
    && chmod +x mosdepth \
    && mv mosdepth /usr/local/bin/mosdepth

RUN command -v mosdepth &>/dev/null || exit 1

COPY [ "requirements.txt", "/opt" ]

# Install Python dependencies
RUN pip install --no-cache-dir -r /opt/requirements.txt

RUN command -v bam2plot &>/dev/null || exit 1

ENTRYPOINT [ "/bin/bash" ]
