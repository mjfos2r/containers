FROM debian:bookworm-slim

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="fasttree" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/fasttree" \
    org.opencontainers.image.description="Debian bookworm-slim container for fasttree" \
    org.opencontainers.image.version="0.0.1"

# Install deps
RUN apt-get update && apt-get install -y \
    ca-certificates \
    curl \
    make \
    fasttree=2.1.11-2 \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Verify installation
RUN for cmd in fasttree; do \
    command -v cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]