FROM debian:bookworm-slim

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="rapidnj" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/rapidnj" \
    org.opencontainers.image.description="Debian bookworm-slim container for rapidnj" \
    org.opencontainers.image.version="2.3.3"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    make \
    # Add additional dependencies here
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Setup directories
RUN mkdir -p /usr/local/src/rapidnj-2.3.3

# Download and install tool
RUN curl -LsSf https://github.com/somme89/rapidNJ/archive/refs/tags/latest.tar.gz \
    | tar xzv -C /usr/local/src/rapidnj-2.3.3 --strip-components=1 --no-same-owner \
    && cd /usr/local/src/rapidnj-2.3.3 \
    && make \
    && cp -R bin/* /usr/local/bin/

# Verify installation
RUN for cmd in rapidNJ; do \
    command -v cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]