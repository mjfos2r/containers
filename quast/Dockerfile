FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/quast" \
    org.opencontainers.image.description="container for Quast v5.3.1" \
    org.opencontainers.image.version="5.3.1"

# swiped from staphb's container.
RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    zlib1g-dev \
    pkg-config \
    libfreetype6-dev \
    libpng-dev \
    wget \
    g++ \
    make \
    perl \
    python3 \
    python3-setuptools \
    python3-matplotlib \
    python3-simplejson \
    python3-joblib \
    locales \
    python3-pip && \
    locale-gen en_US.UTF-8 && \
    apt-get clean && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/*

# Okay so apparently python is hardcoded into the shebang of every single component script. 
# alias me up
RUN echo 'alias python="python3"' >> ~/.bashrc
SHELL [ "/bin/bash", "-c" ] 
RUN mkdir /usr/local/src/quast
##/tmp/bwa /tmp/minimap2
# forgot to bump version so we're re-building and pushing a 5.3.0 image.
## Since busco v3/odb9 dbs are no longer available, use my personal fork instead with the fresh v5 dbs.  
#RUN wget -qO- https://github.com/ablab/quast/releases/download/quast_5.3.0/quast-5.3.0.tar.gz |\
RUN curl -L https://github.com/mjfos2r/quast/raw/refs/heads/master/quast-5.3.1.tar.gz |\
    tar -xz -C /usr/local/src/quast --strip-components=1 --no-same-owner \
    && cd /usr/local/src/quast \
    && ./setup.py install_full

# Verify installation
RUN command -v quast.py  &>/dev/null || exit 1

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash", "--login" ]
