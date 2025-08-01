FROM debian:bookworm-slim
# uv or bust
COPY --from=ghcr.io/astral-sh/uv:0.7.22 /uv /uvx /bin/

# Install dependencies
RUN apt-get update && apt-get install -y \
    bash \
    build-essential \
    git \
    make \
    libomp-dev \
    wget \
    rsync \
    curl \
    ca-certificates \
    perl \
    zlib1g-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# setup UV stuff
RUN uv venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN uv pip install biopython

# Create temp dirs and db dirs
RUN mkdir -p /tmp/kraken2 /tmp/KrakenTools /usr/local/KronaTools /kraken2_dbs /krona_db

# Download kraken2
#RUN wget -qO - https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz | \
#    tar -xz -C /tmp/kraken2 --strip-components=1 --no-same-owner
RUN curl -LsSf https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.6.tar.gz | \
    tar -xzvf - -C /tmp/kraken2 --strip-components=1 --no-same-owner

# Download KronaTools
RUN wget -qO - https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar | \
    tar -x -C /usr/local/KronaTools --strip-components=1 --no-same-owner

# Download KrakenTools
RUN wget -qO - https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.tar.gz | \
    tar -xz -C /tmp/KrakenTools --strip-components=1 --no-same-owner

# Install kraken2
RUN cd /tmp/kraken2 && \
    ./install_kraken2.sh /usr/local/bin/

# Install KrakenTools
RUN cd /tmp/KrakenTools && \
    for script in *.py; do \
        sed -i '1s|^#!.*|#!/usr/bin/env python3|' $script && \
        cp $script /usr/local/bin/ && \
        chmod +x /usr/local/bin/$script; \
        done

# Install KronaTools
RUN cd /usr/local/KronaTools && \
    ./install.pl --taxonomy /krona_db

#RUN cd /tmp/KronaTools && \
#    for script in scripts/*.pl; do \
#        name=$(basename $script .pl); \
#        cp $script /usr/local/bin/kt$name && \
#        chmod +x /usr/local/bin/kt$name; \
#    done && \
#    cp lib/KronaTools.pm /etc/perl/ && \
#    chmod 644 /etc/perl/KronaTools.pm && \
#    cd /

# Test that our tools installed correctly and are in our path!
RUN for CMD in ktImportText kreport2krona.py kraken2 k2; do \
    command -v $CMD &>/dev/null || exit 1; \
    done

# Clean up
RUN rm -rf /tmp/kraken2 /tmp/KrakenTools

# Set environment variables
ENV KRAKEN2_DB_PATH=/kraken2_dbs
ENV PATH=$PATH:/usr/local/bin/kraken2/

ENTRYPOINT [ "/bin/bash" ]
