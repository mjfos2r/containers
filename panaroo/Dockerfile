FROM python:3.10-slim-bookworm
COPY --from=ghcr.io/astral-sh/uv:0.9.15 /uv /uvx /bin/

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="panaroo" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/panaroo" \
    org.opencontainers.image.description="UV managed python container for panaroo" \
    org.opencontainers.image.version="0.0.1"

# installed tools are available immediately
ENV UV_TOOL_BIN_DIR=/usr/local/bin
# enable bytecode compilation
ENV UV_COMPILE_BYTECODE=1

# Install system deps
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    ca-certificates \
    zlib1g-dev \
    libomp-dev \
    make \
    # add more deps here
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# Install cd-hit
RUN mkdir -p /tmp/cd-hit \
    && curl -lsSf https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz \
    | tar xzvf - -C /tmp/cd-hit --strip-components=1 --no-same-owner \
    && cd /tmp/cd-hit \
    && make install \
    && cd - \
    && rm -rf /tmp/cd-hit

# Install mafft
RUN mkdir -p /tmp/mafft \
    && curl -lsSf https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz \
    | tar xzvf - -C /tmp/mafft --strip-components=1 --no-same-owner \
    && cd /tmp/mafft \
    && make clean \
    && make \
    && make install \
    && cd - \
    && rm -rf /tmp/cd-hit

# Download mash binary
RUN curl -lsSf https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar \
    | tar -xv --strip-components=1 --no-same-owner -C /usr/local/bin \
    && chmod +x /usr/local/bin/mash \
    && rm /usr/local/bin/LICENSE.txt

# Download Prank binary
RUN mkdir -p /tmp/prank \
    && curl -lsSf https://github.com/ariloytynoja/prank-msa/blob/master/binaries/prank.linux64.251117.tgz \
    | tar xzvf - -C /tmp/prank --strip-components=1 --no-same-owner \
    && cp /tmp/prank/bin/* /usr/local/bin/

# install panaroo from pip
RUN uv tool install git+https://github.com/gtonkinhill/panaroo

# Verify installation
RUN for cmd in cdhit mafft mash prank panaroo; do \
    command -v cmd &>/dev/null || exit 1; \
    done

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]
