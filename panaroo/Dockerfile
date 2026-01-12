FROM mambaorg/micromamba:1.5.10

# OCI formatted labels
LABEL \
    org.opencontainers.image.title="panaroo" \
    org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.url="https://github.com/mjfos2r/containers/panaroo" \
    org.opencontainers.image.description="Micromamba container for panaroo with MSA tools" \
    org.opencontainers.image.version="0.1.0"

USER root

# Install minimal system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install panaroo and all deps from bioconda as root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.10 \
    panaroo \
    cd-hit \
    mafft \
    mash \
    prank \
    clustalo \
    prokka \
    && micromamba clean --all --yes

WORKDIR /data

ENTRYPOINT ["/bin/bash", "-c"]
