FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/nextdenovo" \
    org.opencontainers.image.description="Container with Nextdenovo and Nextpolish" \
    org.opencontainers.image.version="2.5.2"

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    python3 \
    python3-pip \
    liblzma-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /data /opt/nextdenovo /opt/nextpolish

WORKDIR /opt

RUN pip install paralleltask --break-system-packages

# Download and install NextDenovo
RUN wget -qO- https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz |\
    tar -xvz --strip-components=1 --no-same-owner -C /opt/nextdenovo

# Download and install NextPolish
RUN wget -qO- https://github.com/Nextomics/NextPolish/releases/download/v1.4.1/NextPolish.tgz |\
    tar -xvz --strip-components=1 --no-same-owner -C /opt/nextpolish

# Verify installation
RUN for cmd in nextdenovo nextpolish; do \
    command -v /opt/cmd cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]