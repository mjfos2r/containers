FROM python:3.12-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/username/repository" \
    org.opencontainers.image.description="Python3.12-slim container with ONT Medaka and Samtools for assembly polishing" \
    org.opencontainers.image.version="1.0.0"

RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    make \
    wget \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libdeflate-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

RUN mkdir -p /usr/local/src/htslib /usr/local/src/samtools /tmp/


# Download and compile htslib
RUN wget -qO- https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 \
| tar -xj -C /usr/local/src/htslib --strip-components=1 --no-same-owner \
&& cd /usr/local/src/htslib \
&& ./configure --prefix=/usr/local/ \
&& make \
&& make install

# Download and install samtools
RUN wget -qO- https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
&& tar -xj -C /usr/local/src/samtools --strip-components=1 \
&& cd /usr/local/src/samtools \
&& ./configure --prefix=/usr/local \
&& make \
&& make test \
&& make install

WORKDIR /opt

# Download minimap2 and install to path
RUN wget -qO- https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    | tar xj -C /tmp/minimap2 --strip-components=1 --no-same-owner \
    && mv /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && mv /tmp/minimap2/k8 /usr/local/bin/k8 \
    && mv /tmp/minimap2/paftools.js /usr/local/bin/paftools.js \
    && mv /tmp/minimap2/minimap2.1 /usr/local/share/man/minimap2.1 \
    && rm /tmp/minimap2/* \
    && rmdir /tmp/minimap2

# Set working directory
WORKDIR /opt

# Install Python dependencies
RUN pip install --no-cache-dir medaka

# Test that our tools installed correctly and are in our path!
RUN for cmd in samtools minimap2 tabix medaka; do \
    command -v cmd &>/dev/null || exit 1; \
    done

ENTRYPOINT [ "/bin/bash" ]
