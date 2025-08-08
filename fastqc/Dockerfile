FROM debian:bookworm-slim
COPY --from=ghcr.io/astral-sh/uv:0.7.22 /uv /uvx /bin/

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    ca-certificates \
    unzip \
    gnuplot \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    parallel \
    man-db \
    perl \
    default-jre \
    xxd \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

# install fastqc and symlink to path.
RUN curl -LsSf https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip \
    -o fastqc.zip \
    && unzip fastqc.zip -d /usr/local/src/ \
    && ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc \
    && rm fastqc.zip

# Download and install samtools
RUN mkdir -p /tmp/samtools \
    && curl -LsSf https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    -o /tmp/samtools/samtools.tar.bz2 \
    && tar -xjf /tmp/samtools/samtools.tar.bz2 -C /tmp/samtools --strip-components=1 --no-same-owner \
    && cd /tmp/samtools \
    && ./configure --prefix=/usr/local \
    && make \
    && make test \
    && make install \
    && rm -rf /tmp/samtools

# Download SeqKit to PATH
RUN curl -LsSf https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz \
    | tar -xz -C /usr/local/bin --no-same-owner

# make our venv.
RUN uv venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:$PATH"
# Now install multiqc using uv :)
RUN uv pip install multiqc

# Verify installation
RUN for cmd in fastqc multiqc; do \
    command -v cmd &>/dev/null || exit 1; \
    done

# Set entrypoint
ENTRYPOINT [ "/bin/bash" ]
