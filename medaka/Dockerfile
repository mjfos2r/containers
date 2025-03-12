FROM python:3.12-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/username/repository" \
    org.opencontainers.image.description="Python3.12-slim container with ONT Medaka and Samtools for assembly polishing" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="mfoster11@mgh.harvard.edu"

    RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    make \
    wget \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libdeflate-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Download and install samtools
RUN mkdir -p /tmp/samtools \
    && wget -O /tmp/samtools/samtools.tar.bz2 \
    https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    && tar -xjf /tmp/samtools/samtools.tar.bz2 -C /tmp/samtools --strip-components=1 \
    && cd /tmp/samtools \
    && ./configure --prefix=/usr/local \
    && make \
    && make test \
    && make install \
    && rm -rf /tmp/samtools

# Set working directory
WORKDIR /opt

# Install Python dependencies
RUN pip install --no-cache-dir medaka

RUN command -v samtools &>/dev/null || exit 1
RUN command -v medaka &>/dev/null || exit 1

ENTRYPOINT [ "/bin/bash" ]