# De-conda, reproducible long-read toolchain
FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Base build/runtime deps (C/C++ toolchain, zlib, curl, Python, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git cmake ninja-build pkg-config \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4 ca-certificates \
    python3 python3-pip python3-venv \
    wget curl xz-utils tar gzip parallel \
    && rm -rf /var/lib/apt/lists/*

# Put /usr/local/bin first
ENV PATH="/usr/local/bin:${PATH}"
# plassembler DB location kept from your original
ENV PLASSEMBLER_DB="/usr/local/share/plassembler/db"

# ---------- Assemblers & polishing tools (from upstream repos) ----------
# 1) minimap2 (lh3)  2) miniasm (lh3)
RUN git clone --depth 1 https://github.com/lh3/minimap2.git /tmp/minimap2 && \
    make -C /tmp/minimap2 && cp /tmp/minimap2/minimap2 /usr/local/bin/ && \
    git clone --depth 1 https://github.com/lh3/miniasm.git /tmp/miniasm && \
    make -C /tmp/miniasm && cp /tmp/miniasm/miniasm /usr/local/bin/ && \
    rm -rf /tmp/minimap2 /tmp/miniasm
#   Sources: minimap2 build instructions and repo (lh3). :contentReference[oaicite:0]{index=0}
#            miniasm repo. :contentReference[oaicite:1]{index=1}

# 3) wtdbg2 (redbean)
RUN git clone --depth 1 https://github.com/ruanjue/wtdbg2 /tmp/wtdbg2 && \
    make -C /tmp/wtdbg2 && \
    cp /tmp/wtdbg2/wtdbg2 /tmp/wtdbg2/wtpoa-cns /usr/local/bin/ && \
    rm -rf /tmp/wtdbg2
#   Source & build notes. :contentReference[oaicite:2]{index=2}

# 4) Racon (consensus polishing; official repo now under lbcb-sci)
RUN git clone --recursive --depth 1 https://github.com/lbcb-sci/racon.git /tmp/racon && \
    cmake -S /tmp/racon -B /tmp/racon/build -DCMAKE_BUILD_TYPE=Release && \
    cmake --build /tmp/racon/build -j && \
    cp /tmp/racon/build/bin/racon /usr/local/bin/ && \
    rm -rf /tmp/racon
#   Build instructions (CMake). :contentReference[oaicite:3]{index=3}

# 5) Raven assembler
RUN git clone --depth 1 https://github.com/lbcb-sci/raven /tmp/raven && \
    cmake -S /tmp/raven -B /tmp/raven/build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release && \
    cmake --build /tmp/raven/build -j && \
    cp /tmp/raven/build/raven /usr/local/bin/ && \
    rm -rf /tmp/raven
#   Build instructions. :contentReference[oaicite:4]{index=4}

# 6) Flye (Python package from upstream)
RUN pip3 install --no-cache-dir flye
#   Project home & releases. :contentReference[oaicite:5]{index=5}

# 7) Canu (release tarball)
#    v2.3 is current (Dec 17, 2024); repo recommends binary tarballs.
RUN set -e; \
    cd /tmp; \
    curl -L -o canu.tar.xz https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz && \
    tar -xJf canu.tar.xz && mv canu-2.3/bin/* /usr/local/bin/ && \
    rm -rf /tmp/canu*
#   Releases page and binary install guidance. :contentReference[oaicite:6]{index=6}

# 8) LJA (long-read De Bruijn assembler)
RUN git clone --depth 1 https://github.com/AntonBankevich/LJA /tmp/lja && \
    cmake -S /tmp/lja -B /tmp/lja/build -DCMAKE_BUILD_TYPE=Release && \
    cmake --build /tmp/lja/build -j && \
    cp /tmp/lja/build/bin/* /usr/local/bin/ || true && \
    rm -rf /tmp/lja
#   Repo. (LJA build uses CMake; binary name may vary by version.) :contentReference[oaicite:7]{index=7}

# 9) metaMDBG (metagenome assembler; requires minimap2)
RUN git clone --depth 1 https://github.com/GaetanBenoitDev/metaMDBG /tmp/metamdbg && \
    cmake -S /tmp/metamdbg -B /tmp/metamdbg/build -DCMAKE_BUILD_TYPE=Release && \
    cmake --build /tmp/metamdbg/build -j && \
    cp /tmp/metamdbg/build/bin/metaMDBG /usr/local/bin/ && \
    rm -rf /tmp/metamdbg
#   Build from source option described in README. :contentReference[oaicite:8]{index=8}

# 10) MiniPolish (Python; requires minimap2 & racon present)
RUN pip3 install --no-cache-dir git+https://github.com/rrwick/Minipolish.git
#   Install instructions. :contentReference[oaicite:9]{index=9}

# 11) NECAT (binary distribution from repo)
RUN set -e; \
    cd /tmp; \
    curl -L -o NECAT.tar.gz https://github.com/xiaochuanle/NECAT/archive/refs/heads/master.tar.gz && \
    tar -xzf NECAT.tar.gz && \
    cd NECAT-* && \
    make -j || true && \
    # NECAT ships scripts/binaries; copy common entrypoints if present
    find . -maxdepth 2 -type f -perm -111 -name "necat.pl" -o -name "extract_sequences" -o -name "correct" -o -name "assemble" -o -name "bridge" | xargs -I{} cp {} /usr/local/bin/ || true && \
    rm -rf /tmp/NECAT*
#   Repo (binary/script tool; build targets vary). :contentReference[oaicite:10]{index=10}

# 12) NextDenovo (release tarball + make)
RUN set -e; \
    cd /opt && \
    wget -q https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz && \
    tar -xzf NextDenovo.tgz && rm NextDenovo.tgz && \
    make -C NextDenovo && ln -s /opt/NextDenovo/nextDenovo /usr/local/bin/nextDenovo
#   README shows downloading latest NextDenovo.tgz and make. :contentReference[oaicite:11]{index=11}

# 13) NextPolish (release tarball + make) + its Python dependency
RUN pip3 install --no-cache-dir paralleltask && \
    cd /opt && \
    wget -q https://github.com/Nextomics/NextPolish/releases/download/v1.4.1/NextPolish.tgz && \
    tar -xzf NextPolish.tgz && rm NextPolish.tgz && \
    make -C NextPolish && ln -s /opt/NextPolish/nextPolish /usr/local/bin/nextPolish
#   README shows wget NextPolish.tgz, make; needs paralleltask. :contentReference[oaicite:12]{index=12}

# 14) Flye’s friends already installed; proceed with raven/miniasm/racon handled above.

# 15) NextPolish2 is optional; not requested, so omitted.

# 16) Unicycler (Python; your requested 0.5.0 pin)
RUN pip3 install --no-cache-dir "unicycler==0.5.0"
#   Project home & docs. :contentReference[oaicite:13]{index=13}

# 17) plassembler (Python; your requested 1.8.0 pin)
RUN pip3 install --no-cache-dir "plassembler==1.8.0"
#   Project and docs. :contentReference[oaicite:14]{index=14}

# ---------- Autocycler (binary) ----------
# Pull the latest Linux x86_64 musl tarball automatically and place binary on PATH.
RUN set -e; \
    cd /tmp; \
    curl -s https://api.github.com/repos/rrwick/Autocycler/releases/latest \
     | grep -Eo '"browser_download_url": *"[^"]+linux-x86_64-musl[^"]+tar\.gz"' \
     | cut -d'"' -f4 | head -n1 | xargs -I{} curl -L -o autocycler.tar.gz "{}" && \
    tar -xzf autocycler.tar.gz && \
    mv autocycler /usr/local/bin/ && \
    rm -f autocycler.tar.gz
#   Autocycler latest release (as of Aug 2025: v0.5.2) and binary naming. :contentReference[oaicite:15]{index=15}

# ---------- housekeeping ----------
# GNU parallel citation prompt opt-out (per its docs)
RUN touch /root/.parallel && chmod 600 /root/.parallel

# Create and own the plassembler DB dir, then fetch DB
RUN mkdir -p /usr/local/share/plassembler && \
    plassembler download -d "${PLASSEMBLER_DB}"

# Add your entrypoint
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

