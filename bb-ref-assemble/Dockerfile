FROM python:3.13-slim
COPY --from=ghcr.io/astral-sh/uv:0.7.8 /uv /uvx /bin/

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    ca-certificates \
    git \
    nano \
    vim \
    gnuplot \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    parallel \
    man-db \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Setup all of our directories
RUN mkdir -p /opt/uv-bin /tmp/minimap2 /usr/local/src/samtools /usr/local/src/htslib /usr/local/src/racon /usr/local/src/rebaler

# and uv's tool box
ENV UV_TOOL_BIN_DIR="/opt/uv-bin/"
ENV PATH="/opt/uv-bin:$PATH"

# install seqkit
RUN curl -LsS https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz >/tmp/seqkit_linux_amd64.tar.gz \
    && curl -LsS https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz.md5.txt >/tmp/seqkit_linux_amd64.tar.gz.md5.txt \
    && cd /tmp && md5sum -c seqkit_linux_amd64.tar.gz.md5.txt \
    && tar -xz -C /usr/local/bin --no-same-owner -f seqkit_linux_amd64.tar.gz

# install minimap2
RUN curl -LsS https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    | tar xj -C /tmp/minimap2 --strip-components=1 --no-same-owner \
    && mv /tmp/minimap2/minimap2 /usr/local/bin/minimap2 \
    && mv /tmp/minimap2/k8 /usr/local/bin/k8 \
    && mv /tmp/minimap2/paftools.js /usr/local/bin/paftools.js \
    && mv /tmp/minimap2/minimap2.1 /usr/local/share/man/minimap2.1 \
    && rm /tmp/minimap2/* \
    && rmdir /tmp/minimap2

# Get and compile htslib
RUN curl -LsS https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 \
    | tar -xj -C /usr/local/src/htslib --strip-components=1 --no-same-owner \
    && cd /usr/local/src/htslib \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install

# Get samtools
RUN curl -LsS https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    | tar -xj -C /usr/local/src/samtools --strip-components=1 --no-same-owner \
    && cd /usr/local/src/samtools \
    && ./configure --prefix=/usr/local/ \
    && make \
    && make install

# get racon
RUN cd /usr/local/src && git clone --recursive https://github.com/lbcb-sci/racon.git racon \
    && cd racon \
    && mkdir build && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release ..\
    && make install

# get rebaler
RUN cd /usr/local/src && git clone --recursive https://github.com/rrwick/Rebaler.git rebaler \
    && cd rebaler \
    && uv tool install .

ENTRYPOINT [ "/bin/bash" ]