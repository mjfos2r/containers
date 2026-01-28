FROM python:3.11-slim-bookworm
COPY --from=ghcr.io/astral-sh/uv:0.9.15 /uv /uvx /bin/

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="atc_clean_gfa" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/atc_clean_gfa" \
    org.opencontainers.image.description="UV managed python container for atc_clean_gfa" \
    org.opencontainers.image.version="0.0.1"

# installed tools are available immediately
ENV UV_TOOL_BIN_DIR=/usr/local/bin
# enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    ca-certificates \
    wget \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Install Rust and cargo
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:$PATH"

# Install autocycler 0.5.2 via cargo
RUN cargo install autocycler --version 0.5.2

# Install seqkit
RUN wget -qO- https://github.com/shenwei356/seqkit/releases/download/v2.8.0/seqkit_linux_amd64.tar.gz \
    | tar -xzf - -C /usr/local/bin/

# Copy the cleaning script
COPY clean_gfa.py /usr/local/bin/clean_gfa.py
RUN chmod +x /usr/local/bin/clean_gfa.py

# Verify installation
RUN for cmd in autocycler seqkit clean_gfa.py; do \
    command -v $cmd &>/dev/null || exit 1; \
    done

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]