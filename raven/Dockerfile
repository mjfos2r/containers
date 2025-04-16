FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/raven" \
    org.opencontainers.image.description="Container for Raven assembler" \
    org.opencontainers.image.version="1.8.3"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    cmake \
    git \
    python3 \
    python3-pip \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Clone and install Raven
RUN git clone https://github.com/lbcb-sci/raven.git \
    && cd raven \
    && cmake -S ./ -B./build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build \
    && cmake --install ./build

# install python bindings
RUN pip install git+git://github.com/lbcb-sci/raven.git@master

# Verify installation
RUN command -v raven &>/dev/null || exit 1

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]