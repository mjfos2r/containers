FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/necat" \
    org.opencontainers.image.description="Container for NECAT assembler" \
    org.opencontainers.image.version="0.0.1"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    perl \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# Download NECAT to path
RUN wget -qO- https://github.com/xiaochuanle/NECAT/releases/download/v0.0.1_update20200803/necat_20200803_Linux-amd64.tar.gz \
| tar -xz --strip-components=3 --no-same-owner -C /usr/local/bin

# Verify installation
RUN for cmd in config correct assemble bridge; do \
    command -v necat.pl cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]