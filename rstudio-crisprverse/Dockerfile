FROM rocker/rstudio:4.5.1
RUN /rocker_scripts/install_tidyverse.sh
RUN /rocker_scripts/install_tensorflow.sh
RUN apt-get install -y zlib1g-dev

RUN echo "Downloading gcc-6-base" \
    && cd /tmp/ && wget http://archive.ubuntu.com/ubuntu/pool/universe/g/gcc-6/gcc-6-base_6.4.0-17ubuntu1_amd64.deb \
    && echo "Downloading libgfortran3" \ 
    && cd /tmp/ && wget http://archive.ubuntu.com/ubuntu/pool/universe/g/gcc-6/libgfortran3_6.4.0-17ubuntu1_amd64.deb \
    && echo "Installing gcc-6-base" \
    && cd /tmp/ && sudo dpkg -i gcc-6-base_6.4.0-17ubuntu1_amd64.deb \ 
    && echo "Installing libgfortran3" \
    && cd /tmp/ && sudo dpkg -i libgfortran3_6.4.0-17ubuntu1_amd64.deb

RUN R -e "BiocManager::install('crisprVerse')"
RUN R -e "BiocManager::install('crisprBwa')"
RUN R -e "BiocManager::install('BSgenomeForge')"
RUN R -e "BiocManager::install('devtools')"
RUN R -e "BiocManager::install('txdbmaker')"

COPY --chown=rstudio:rstudio rstudio-prefs.json /home/rstudio/.config/rstudio
RUN mkdir /home/rstudio/.cache/R/ExperimentHub
