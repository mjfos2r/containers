FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/quast" \
    org.opencontainers.image.description="container for Quast v5.3.0" \
    org.opencontainers.image.version="5.3.0"

# swiped from staphb's container.
RUN apt-get update && apt-get install -y \
    zlib1g-dev \
    pkg-config \
    libfreetype6-dev \
    libpng-dev \
    wget \
    g++ \
    make \
    perl \
    python \
    python-setuptools \
    python-matplotlib \
    python-simplejson \
    python-joblib \
    locales \
    python-pip && \
    locale-gen en_US.UTF-8 && \
    apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /tmp/quast /tmp/bwa /tmp/minimap2
# Download and install tool
RUN wget -qO- https://github.com/ablab/quast/releases/download/quast_5.3.0/quast-5.3.0.tar.gz |\
    tar -xz - -C /tmp/quast --strip-components=1 --no-same-owner \
    && cd /tmp/quast \
    && python3 setup.py install_full

# Verify installation
RUN command -v quast.py  &>/dev/null || exit 1

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]