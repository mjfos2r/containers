FROM mambaorg/micromamba:latest

# I guess I'll authenticate gurobi each time I instantiate this container.
# I'll figure out how to streamline this for terra if this turns out to be
# worthwhile.
RUN micromamba create -y -n phables python=3.12 -c conda-forge -c bioconda -c anaconda && \
    for tool in phables gurobi flye; do \
        if micromamba install -y -n autocycler -c conda-forge -c bioconda "$tool"; then \
            printf " $tool" >> /tmp/installed_tools.txt; \
        else \
            printf " $tool" >> /tmp/missing_tools.txt; \
        fi \
    done && \
    printf "\n" >> /tmp/installed_tools.txt && printf "\n" >> /tmp/missing_tools.txt && \
    micromamba clean --all --yes

# Add the autocycler environment to PATH.
ENV PATH="/opt/conda/envs/phables/bin:$PATH"

ENTRYPOINT [ "/bin/bash" ]
