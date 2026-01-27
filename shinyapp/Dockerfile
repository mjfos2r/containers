FROM rocker/shiny:4.5.2

# Install deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# install r-packages from CRAN
RUN R -e "install.packages( \
            c('shiny', 'shinyWidgets', 'ape', 'ggplot2', \
            'ggnewscale', 'ggiraph', 'BiocManager', \
            'dplyr', 'tidyr', 'phytools', 'showtext'), \
            repos = 'https://cloud.r-project.org' \
            )"

# install bioconductor packages
RUN R -e "BiocManager::install(c('ggtree', 'ggtreeExtra'), ask=FALSE, update=FALSE)"

WORKDIR /srv/shiny-server/InteractiveBorreliaPhylogeny

COPY app.R .
COPY data ./data

# ensure user: shiny can see and run everything
RUN chown -R shiny:shiny /srv/shiny-server

# expose shiny's port
EXPOSE 3838

# use default shiny-server cmd from base image
CMD ["/usr/bin/shiny-server"]
