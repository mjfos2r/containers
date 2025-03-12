FROM python:3.12-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/username/repository" \
    org.opencontainers.image.description="Python3.12-slim container for Flye longread assembler" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="<mfoster11@mgh.harvard.edu>"

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    g++ \
    make \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/flye

WORKDIR /tmp

RUN pip install --no-cache-dir setuptools

RUN git clone https://github.com/mikolmogorov/Flye.git \
    && cd Flye \
    && python setup.py install --prefix="/usr/local/" \
    && cd .. \
    && rm -rf /tmp/Flye

RUN command -v flye &>/dev/null || exit 1

ENTRYPOINT ["python", "/opt/main.py"]
