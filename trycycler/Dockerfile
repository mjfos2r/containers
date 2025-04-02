FROM python:3.11-slim

# Metadata labels using OCI standard format
LABEL org.opencontainers.image.authors="<mfoster11@mgh.harvard.edu>" \
    org.opencontainers.image.source="https://github.com/mjfos2r/containers/trycycler" \
    org.opencontainers.image.description="python3.11-slim container with Trycycler assembler" \
    org.opencontainers.image.version="1.0.0"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    r-base \
    r-base-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /tmp/muscle /tmp/minimap2

# Set working directory
WORKDIR /tmp

# clone and install miniasm (Trycycler dep 1)
RUN git clone https://github.com/lh3/miniasm.git  \
    && cd miniasm \
    && make

# Get minimap2 and install to path
RUN wget -qO- https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    | tar -xj -C /tmp/minimap2 --strip-components=1 --no-same-owner \
    && mv /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && mv /tmp/minimap2/k8 /usr/local/bin/k8 \
    && mv /tmp/minimap2/paftools.js /usr/local/bin/paftools.js \
    && mv /tmp/minimap2/minimap2.1 /usr/local/share/man/minimap2.1 \
    && rm /tmp/minimap2/* \
    && rmdir /tmp/minimap2

# Download mash binary to PATH and cleanup.
RUN wget -qO- https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar \
    | tar -xv --strip-components=1 --no-same-owner -C /usr/local/bin\
    && chmod +x /usr/local/bin/mash \
    && rm /usr/local/bin/LICENSE.txt

# download muscle_v3, compile, and move binary to path.
RUN wget -qO- https://drive5.com/muscle/muscle_src_3.8.1551.tar.gz \
    | tar -xz --no-same-owner -C /tmp/muscle \
    && cd /tmp/muscle \
    && make \
    && cd ../ \
    && mv /tmp/muscle/muscle /usr/local/bin/muscle \
    && rm -rf /tmp/muscle

# Install required R packages for trycycler cluster step.
RUN R -e "install.packages(c('ape', 'phangorn'), repos='https://cloud.r-project.org/')"
# Verify R installation
RUN R -e "if (!all(c('ape', 'phangorn') %in% installed.packages()[,'Package'])) { quit(status=1) }"

WORKDIR /opt

# Install trycycler via pip
RUN git clone https://github.com/rrwick/Trycycler.git \
    && pip install ./Trycycler \
    && cp Trycycler/scripts/canu_trim.py /usr/local/bin/

# Verify the rest of our installations
RUN for cmd in minimap2 mash muscle trycycler; do \
    command -v cmd &>/dev/null || exit 1; \
    done

RUN command -v python3 /usr/local/bin/canu_trim.py &>/dev/null || exit 1;

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]