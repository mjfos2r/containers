FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="mfoster11@mgh.harvard.edu" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/mummer4" \
    org.opencontainers.image.description="container for mummer4" \
    org.opencontainers.image.version="1.1.0" \
    maintainer="mfoster11@mgh.harvard.edu"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    perl \
    bash \
    gnuplot \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Setup all of our directories
RUN mkdir -p /data/ref /alignments /output /opt/mummer
# get mummer4 ...pls work...
RUN wget -P /opt/mummer -O /opt/mummer/mummer.tar.gz https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
RUN tar -xzf /opt/mummer/mummer.tar.gz -C /opt/mummer --strip-components=1 --no-same-owner \
    && cd /opt/mummer \
    && ./configure LDFLAGS=-static \
    && make \
    && make install

# Test that our tools installed correctly and are in our path!
RUN for CMD in mummer nucmer promer; do \
    command -v $CMD &>/dev/null || echo "$CMD not found!" && exit 1; \
    done

# Change context to /opt
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]