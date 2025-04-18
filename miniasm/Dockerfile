FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/miniasm" \
    org.opencontainers.image.description="Container for miniasm and minipolish" \
    org.opencontainers.image.version="0.3.0"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    gcc \
    zlib1g-dev \
    python3 \
    python3-pip \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN mkdir /tmp/minimap2

# Get minimap2 and install to path
RUN wget -qO- https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    | tar xj -C /tmp/minimap2 --strip-components=1 --no-same-owner \
    && mv /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && mv /tmp/minimap2/k8 /usr/local/bin/k8 \
    && mv /tmp/minimap2/paftools.js /usr/local/bin/paftools.js \
    && mv /tmp/minimap2/minimap2.1 /usr/local/share/man/minimap2.1 \
    && rm /tmp/minimap2/* \
    && rmdir /tmp/minimap2

# clone and install miniasm
RUN git clone https://github.com/lh3/miniasm.git  \
    && cd miniasm \
    && make

# clone and install minipolish
RUN git clone https://github.com/rrwick/Minipolish.git \
    && pip3 install ./Minipolish --break-system-packages \
    && mv Minipolish/miniasm_and_minipolish.sh /usr/local/bin/miniasm_and_minipolish.sh

# install any2fasta
RUN wget https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta -o /usr/local/bin \
    && chmod +x /usr/local/bin/any2fasta

# Verify installation
RUN for cmd in minimap2 k8 paftools.js miniasm minipolish miniasm_and_minipolish.sh any2fasta; do \
    command -v cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]