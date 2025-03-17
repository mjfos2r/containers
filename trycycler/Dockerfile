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


# Set working directory
WORKDIR /opt

# clone and install miniasm (Trycycler dep 1)
RUN git clone https://github.com/lh3/miniasm.git  \
    && cd miniasm \
    && make

# Get minimap2 and install to path
RUN wget -qO- https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    | tar xj -C /tmp/minimap2 --strip-components=1 --no-same-owner \
    && mv /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && mv /tmp/minimap2/k8 /usr/local/bin/k8 \
    && mv /tmp/minimap2/paftools.js /usr/local/bin/paftools.js \
    && mv /tmp/minimap2/minimap2.1 /usr/local/share/man/minimap2.1 \
    && rm /tmp/minimap2/* \
    && rmdir /tmp/minimap2

# Download mash binary to PATH and cleanup.
RUN wget -qO- https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar \
    | tar -xvf -C /usr/local/bin --strip-components=1 \
    && chmod +x /usr/local/bin/mash \
    && rm /usr/local/bin/LICENSE.txt

# download muscle_v3 to PATH
RUN wget -qO- https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz \
    | tar -xzf -C /usr/local/bin \
    && mv /usr/local/bin/muscle3.8.31_i86linux64 /usr/local/bin/muscle \
    && chmod +x /usr/local/bin/muscle

# Install required R packages for trycycler cluster step.
RUN R -e "install.packages(c('ape', 'phangorn'), repos='https://cloud.r-project.org/')"

# Install trycycler via pip
RUN git clone https://github.com/rrwick/Trycycler.git \
    && pip install ./Trycycler

# Verify R installation
RUN R -e "if (!all(c('ape', 'phangorn') %in% installed.packages()[,'Package'])) { quit(status=1) }"

# Verify the rest of our installations
RUN for cmd in minimap2 mash muscle trycycler; do \
    command -v cmd &>/dev/null || exit 1; \
    done

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]