FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="TBD" \
    org.opencontainers.image.description="Bookworm-slim container with Chopper-0.9.1" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="<mfoster11@mgh.harvard.edu>"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    ca-certificates \
    zlib1g-dev \
    pigz \
    unzip \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# get the chopper
RUN wget -q https://github.com/wdecoster/chopper/releases/download/v0.9.1/chopper-linux.zip -O /tmp/chopper.zip \
    && unzip /tmp/chopper.zip -d /usr/local/bin \
    && rm /tmp/chopper.zip \
    && chmod +x /usr/local/bin/chopper

# Download and Install SeqKit
RUN wget -qO- https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz \
| tar -xz -C /usr/local/bin --no-same-owner \
    && chmod +x /usr/local/bin/seqkit

# Verify installation
RUN command -v chopper &>/dev/null || exit 1
RUN command -v seqkit &>/dev/null || exit 1

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]
