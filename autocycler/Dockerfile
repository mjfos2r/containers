FROM mambaorg/micromamba:latest

# This container was originally written by Ryan Wick and can be found at:
# https://github.com/rrwick/Autocycler/tree/main/pipelines/Dockerfile_by_Ryan_Wick
# It has been adapted by Michael Foster

# Create an "autocycler" environment with long-read assemblers. (added plassembler 1.8.0 and specified unicycler version 0.5.0)
RUN micromamba create -y -n autocycler python=3.12 -c conda-forge -c bioconda curl \
    && for tool in canu flye lja metamdbg miniasm minimap2 minipolish necat nextdenovo nextpolish racon raven-assembler wtdbg plassembler unicycler; do \
        if micromamba install -y -n autocycler -c conda-forge -c bioconda "$tool"; then \
            printf " $tool" >>/tmp/installed_tools.txt; \
        else \
            printf " $tool" >>/tmp/missing_tools.txt; \
        fi \
    done \
    && printf "\n" >> /tmp/installed_tools.txt && printf "\n" >> /tmp/missing_tools.txt \
    && micromamba clean --all --yes

# Add the autocycler environment to PATH.
ENV PATH="/opt/conda/envs/autocycler/bin:$PATH"
# also set the global for plassembler_db
ENV PLASSEMBLER_DB="/usr/local/share/plassembler/db"

# Download and install Autocycler (binary and scripts) from the latest release.
#RUN curl -L -o autocycler.tar.gz $(curl -s https://api.github.com/repos/rrwick/Autocycler/releases/latest | grep "browser_download_url.*linux-x86_64-musl.*tar.gz" | cut -d '"' -f 4) && \
RUN curl -L -o autocycler.tar.gz https://github.com/rrwick/Autocycler/releases/download/v0.5.1/autocycler-linux-x86_64-musl-v0.5.1.tar.gz && \
    tar -xzf autocycler.tar.gz && \
    mv autocycler /opt/conda/envs/autocycler/bin/ && \
    rm autocycler.tar.gz

# change to root and install parallel via apt for the assembly script and set up some other things for later in the build.
USER root
RUN apt-get update && apt-get install -y parallel \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# also set up the dir for our plassembler_db while we're root.
RUN mkdir -p /usr/local/share/plassembler \
    && chown -R $MAMBA_USER:$MAMBA_USER /usr/local/share/plassembler

RUN install -d -m 1777 /usr/local/share/plassembler

# also add our entrypoint script to root
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

# since unicycler vanishes whenever plassembler goes looking, symlink it to global path.
RUN ln -sf /opt/conda/envs/autocycler/bin/unicycler   /usr/local/bin/unicycler   && \
    ln -sf /opt/conda/envs/autocycler/bin/autocycler  /usr/local/bin/autocycler  && \
    ln -sf /opt/conda/envs/autocycler/bin/plassembler /usr/local/bin/plassembler

# Ensure directory traversal for arbitrary UIDs
RUN chmod 755 /opt /opt/conda /opt/conda/envs /opt/conda/envs/autocycler /opt/conda/envs/autocycler/bin \
    && chmod a+rx /opt/conda/envs/autocycler/bin/unicycler /opt/conda/envs/autocycler/bin/python*
# Optional: fail build if unicycler is missing/broken
RUN micromamba run -n autocycler unicycler --version
# add plassembler download and env var setting.
RUN plassembler download -d "$PLASSEMBLER_DB"

# Install Seqkit
RUN curl -LsSf https://github.com/shenwei356/seqkit/releases/download/v2.10.1/seqkit_linux_amd64.tar.gz \
    | tar -xz -C /usr/local/bin --no-same-owner

ENTRYPOINT [ "/entrypoint.sh" ]