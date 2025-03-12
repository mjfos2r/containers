FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="TBD" \
    org.opencontainers.image.description="Basic slim container running on bookworm-slim for general tasks" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="<mfoster11@mgh.harvard.edu>"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    pigz \
    && rm -rf /var/lib/apt/lists/*

# Setup directories
RUN mkdir -p /data /output

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]