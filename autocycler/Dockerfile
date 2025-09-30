FROM mambaorg/micromamba:latest

# This container was originally written by Ryan Wick and can be found at:
# https://github.com/rrwick/Autocycler/tree/main/pipelines/Dockerfile_by_Ryan_Wick
# It has been adapted by Michael Foster

# Create an "autocycler" environment with long-read assemblers. (added plassembler 1.8.0 and specified unicycler version 0.5.0)
RUN micromamba create -y -n autocycler python=3.12 -c conda-forge -c bioconda curl \
    && for tool in canu flye lja metamdbg miniasm minimap2 minipolish necat nextdenovo nextpolish racon raven-assembler wtdbg "plassembler==1.8.1" "unicycler==0.5.1"; do \
        if micromamba install -y -n autocycler -c conda-forge -c bioconda "$tool"; then \
            printf " $tool" >>/tmp/installed_tools.txt; \
        else \
            printf " $tool" >>/tmp/missing_tools.txt; \
        fi \
    done \
    && printf "\n" >> /tmp/installed_tools.txt && printf "\n" >> /tmp/missing_tools.txt \
    && micromamba clean --all --yes

# Put env bin on PATH for *all* users, including root
ENV PATH="/opt/conda/envs/autocycler/bin:$PATH"
ENV PLASSEMBLER_DB="/usr/local/share/plassembler/db"

# Autocycler binary
RUN curl -L -o autocycler.tar.gz https://github.com/rrwick/Autocycler/releases/download/v0.5.1/autocycler-linux-x86_64-musl-v0.5.1.tar.gz \
    && tar -xzf autocycler.tar.gz \
    && mv autocycler /opt/conda/envs/autocycler/bin/ \
    && rm autocycler.tar.gz

# change to root and install parallel via apt for the assembly script and set up some other things for later in the build.
USER root
RUN apt-get update && apt-get install -y parallel \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /usr/local/share/plassembler

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# Make sure key env binaries are visible globally (root, arbitrary UID)
RUN ln -sf /opt/conda/envs/autocycler/bin/unicycler   /usr/local/bin/unicycler   \
 && ln -sf /opt/conda/envs/autocycler/bin/autocycler  /usr/local/bin/autocycler  \
 && ln -sf /opt/conda/envs/autocycler/bin/plassembler /usr/local/bin/plassembler

# Sanity checks at build time (fail if missing)
RUN unicycler --version \
 && autocycler --version \
 && plassembler --version

# Download DB
RUN plassembler download -d "$PLASSEMBLER_DB"

# Seqkit
RUN curl -LsSf https://github.com/shenwei356/seqkit/releases/download/v2.10.1/seqkit_linux_amd64.tar.gz \
    | tar -xz -C /usr/local/bin --no-same-owner

ENTRYPOINT [ "/entrypoint.sh" ]