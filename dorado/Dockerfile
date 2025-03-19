FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/dorado" \
    org.opencontainers.image.description="Container for ONT Dorado" \
    org.opencontainers.image.version="0.9.1" \

# Install deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    git \
    ca-certificates \
    build-essential \
    libhdf5-dev \
    libssl-dev \
    libzstd-dev \
    autoconf \
    automake

# install specific versions and update as per ONT's documentation: DEV.md
RUN apt install gcc-9 g++-9 --no-install-recommends
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-9
# purge the apt cmake and install the specific version to work with dorado as per ONT
RUN sudo apt purge cmake

# Setup directories
RUN mkdir -p /data /output /opt/tool

# Download and install tool
RUN wget -P /opt/tool -O /opt/tool/tool.tar.gz https://example.com/tool/download/url
RUN tar -xzf /opt/tool/tool.tar.gz -C /opt/tool --strip-components=1 --no-same-owner \
    && cd /opt/tool \
    && ./configure --your-configure-options \
    && make \
    && make install

# Verify installation
RUN for cmd in tool1 tool2 tool3; do \
    command -v cmd &>/dev/null || exit 1; \
done

# Set working directory
WORKDIR /opt

COPY [ "entrypoint.sh", "script.sh", "/opt/" ]

# Make scripts executable
RUN chmod +x /opt/entrypoint.sh /opt/script.sh

# Set entrypoint
ENTRYPOINT [ "/bin/bash", "entrypoint.sh" ]