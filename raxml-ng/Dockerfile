FROM debian:bookworm-slim

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="raxmlhpc" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/raxml-ng" \
    org.opencontainers.image.description="Debian bookworm-slim container for raxml-ng v1.2.2" \
    org.opencontainers.image.version="1.2.2"

# Install deps
RUN apt-get update && apt-get install -y \
    ca-certificates \
    curl \
    unzip \
    libomp-dev \
    # Add additional dependencies here
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/raxml
# Download and install tool
RUN curl -LsSf https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip \
    -o /tmp/raxml/raxml-ng.zip \
    && cd /tmp/raxml \
    && unzip raxml-ng.zip \
    && mv raxml-ng /usr/local/bin

# test RAxML
RUN raxml-ng --msa /tmp/raxml/data/testDNA.fa --model GTR+G || exit 1
RUN raxml-ng --all --msa /tmp/raxml/data/testAA.fa --model LG+G8+F --tree pars{10} --bs-trees 200 || exit 1
RUN raxml-ng --evaluate --msa /tmp/raxml/data/testAA.fa --model partitions.txt --tree test.tree --brlen scaled || exit 1

#RUN rm -rf /tmp/raxml

# Verify installation
RUN for cmd in raxml-ng; do \
    command -v cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]