FROM python:3.12-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/username/repository" \
    org.opencontainers.image.description="Python3.12-slim container with ONT Medaka and Samtools for assembly polishing" \
    org.opencontainers.image.version="2.0.1"

RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    make \
    wget \
    curl \
    ca-certificates \
    git \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libdeflate-dev

# prep the repo to install git-lfs
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
RUN apt-get update && apt-get -y \
    git-lfs \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

RUN mkdir /app/medaka
# Set working directory for medaka installation
WORKDIR /app
# install medaka from source using git-lfs as configured above
# htslib, samtools, and minimap2 are also built from source in this process.
RUN git clone https://github.com/nanoporetech/medaka.git \
    && cd medaka \
    && make install

WORKDIR /opt

RUN cat > entrypoint.sh << EOF
#!/bin/bash
source /app/medaka/venv/bin/activate
exec "\$@"
EOF

RUN chmod +x entrypoint.sh

# Test that our tools installed correctly and are in our path!
RUN source /app/medaka/venv/bin/activate \
    && for cmd in samtools minimap2 tabix medaka; do \
    command -v cmd &>/dev/null || exit 1; \
    done

ENTRYPOINT [ "/opt/entrypoint.sh" ]
CMD [ "medaka", "--version" ]
