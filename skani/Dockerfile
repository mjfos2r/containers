FROM debian:bookworm-slim
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/username/repository" \
    org.opencontainers.image.description="container for skani" \
    org.opencontainers.image.version="1.0.0" \
    maintainer="<mfoster11@mgh.harvard.edu>"

WORKDIR /opt

# Install deps
RUN apt-get update && apt-get install -y \
    wget \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN wget -p /usr/local/bin https://github.com/bluenote-1577/skani/releases/download/latest/skani \
    && chmod +x /usr/local/bin/skani

RUN command -v skani &>/dev/null || exit 1

ENTRYPOINT [ "/bin/bash" ]