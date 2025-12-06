FROM debian:bookworm-slim

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="bandageng" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/bandageng" \
    org.opencontainers.image.description="Debian bookworm-slim container for BandageNG" \
    org.opencontainers.image.version="0.0.1"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    qt6-base-dev \
    qt6-svg-dev \
    libxkbcommon-dev \
    # Add additional dependencies here
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Setup directories
RUN mkdir -p /tmp/bandageng

# install cmake 4.2.0
RUN curl -LsSf https://github.com/Kitware/CMake/releases/download/v4.2.0/cmake-4.2.0-linux-x86_64.sh > /tmp/cmake-install.sh \
    && chmod +x /tmp/cmake-install.sh \
    && /tmp/cmake-install.sh --skip-license --prefix=/usr/local \
    && rm /tmp/cmake-install.sh

# Download and install tool
RUN curl -LsSf https://github.com/asl/BandageNG/archive/refs/tags/v2025.12.1.tar.gz \
    | tar xzv -C /tmp/bandageng --strip-components=1 --no-same-owner \
    && cd /tmp/bandageng \
    && mkdir -p build \
    && cd build \
    && cmake .. \
    && make \
    && make install

# Verify installation
RUN for cmd in BandageNG; do \
    command -v cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]