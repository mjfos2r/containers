FROM python:3.11-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers" \
    org.opencontainers.image.description="Container for telomere resolution in Borrelia burgdorferi fastq files" \
    org.opencontainers.image.version="0.2.0" \
    maintainer="mfoster11@mgh.harvard.edu"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Download and Install SeqKit
RUN wget -qO- https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz \
| tar -xz -C /usr/local/bin --no-same-owner \
    && chmod +x /usr/local/bin/seqkit

# Install Python dependencies
RUN pip install --no-cache-dir -r biopython

COPY [ "clip_telomeric_reads.py", "filter_and_merge.py", "/opt"]
# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]