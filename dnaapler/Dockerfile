FROM debian:bookworm-slim
COPY --from=ghcr.io/astral-sh/uv:0.7.22 /uv /uvx /bin/

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    make \
    libomp-dev \
    wget \
    rsync \
    curl \
    ca-certificates \
    perl \
    pigz \
    zlib1g-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# make our venv.
ENV PATH="/root/.local/bin/:$PATH"

# Install MMseqs2 v13.45111
RUN curl -LsSf https://github.com/soedinglab/MMseqs2/releases/download/13-45111/mmseqs-linux-avx2.tar.gz \
    | tar -xzvf - -C /usr/local/src/
# add it to path
ENV PATH="$PATH:/usr/local/src/mmseqs/bin/"
# install dnaapler
RUN uv tool install dnaapler
# make sure everything is installed
RUN for cmd in mmseqs dnaapler; do \
        command -v "$cmd" &>/dev/null || exit 1; \
    done

