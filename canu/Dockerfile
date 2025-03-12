FROM openjdk:8

LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/canu" \
    org.opencontainers.image.description="openjdk8 container for canu assembler 2.3" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="<mfoster11@mgh.harvard.edu>"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    gnuplot \
    xz-utils \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Setup directories
RUN mkdir -p /data

WORKDIR /tmp
# Download and install tool
RUN wget https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz \
    && echo "5f5e537346f21e91393b0e5447f45bb3 canu-2.3.Linux-amd64.tar.xz" | md5sum -c - || exit 1 \
    && tar -xJf /tmp/canu-2.3.Linux-amd64.tar.xz -C /usr/local/ --strip-components=1 --no-same-owner \
    && rm -rf /tmp/canu-2.3*

RUN command -v canu &>/dev/null || exit 1

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]