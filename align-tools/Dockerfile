FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="mfoster11@mgh.harvard.edu" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/align" \
    org.opencontainers.image.description="Container for Minimap2 with samtools, htslib, and seqkit for use in pipelines." \
    org.opencontainers.image.version="0.1.1" \
    maintainer="mfoster11@mgh.harvard.edu"
# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    ca-certificates \
    gnuplot \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    parallel \
    man-db \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Setup all of our directories
RUN mkdir -p /data /opt /tmp/minimap2 /usr/local/src/samtools /usr/local/src/htslib

# Get minimap2 and install to path
RUN wget -qO- https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    | tar xj -C /tmp/minimap2 --strip-components=1 --no-same-owner \
    && mv /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && mv /tmp/minimap2/k8 /usr/local/bin/k8 \
    && mv /tmp/minimap2/paftools.js /usr/local/bin/paftools.js \
    && mv /tmp/minimap2/minimap2.1 /usr/local/share/man/minimap2.1 \
    && rm /tmp/minimap2/* \
    && rmdir /tmp/minimap2

# Get and compile htslib
RUN wget -qO- https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 \
    | tar -xj -C /usr/local/src/htslib --strip-components=1 --no-same-owner \
    && cd /usr/local/src/htslib \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install

# Get and compile samtools
RUN wget -qO- https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    | tar -xj -C /usr/local/src/samtools --strip-components=1 --no-same-owner \
    && cd /usr/local/src/samtools \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install

# Download and Install SeqKit
RUN wget -qO- https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz \
    | tar -xz -C /usr/local/bin --no-same-owner

# Test that our tools installed correctly and are in our path!
RUN for CMD in samtools minimap2 seqkit; do \
    command -v $CMD &>/dev/null || exit 1; \
    done

# Change context to /opt
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]
CMD [ "minimap2", "--help" ]