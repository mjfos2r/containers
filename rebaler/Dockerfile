FROM debian:bookworm-slim
COPY --from=ghcr.io/astral-sh/uv:0.7.22 /uv /uvx /bin/

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    curl \
    ca-certificates \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN mkdir /tmp/minimap2 /tmp/racon

# Get minimap2 and install to path
RUN curl -LsSf https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    | tar xj -C /tmp/minimap2 --strip-components=1 --no-same-owner \
    && mv /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && mv /tmp/minimap2/k8 /usr/local/bin/k8 \
    && mv /tmp/minimap2/paftools.js /usr/local/bin/paftools.js \
    && mv /tmp/minimap2/minimap2.1 /usr/local/share/man/minimap2.1 \
    && rm /tmp/minimap2/* \
    && rmdir /tmp/minimap2

# make our venv.
RUN uv venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# we need to install racon now.
RUN git clone --recursive https://github.com/lbcb-sci/racon.git /tmp/racon
RUN cd /tmp/racon \
    && mkdir build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    && make \
    && mv bin/racon /usr/local/bin/racon \
    && cd /opt \
    && rm -rf /tmp/racon

# now we install rebaler (my fork.)
RUN uv pip install git+https://github.com/mjfos2r/rebaler.git


# Verify installation
RUN for cmd in racon minimap2 rebaler; do \
    command -v cmd &>/dev/null || exit 1; \
    done

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]
