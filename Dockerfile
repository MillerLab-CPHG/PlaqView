# Install R version 4.0.5
FROM r-base:4.0.5

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
    libxml2-dev \
    libnlopt-dev \
    libudunits2-dev \
    libgeos-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    git         \
    libgdal-dev \
    python3-venv
    
# Install Shiny server
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt) && \
    wget "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb
    
    
#### ESSENTIAL R PACKAGES ####
## FROM CRAN 
RUN install2.r --error --skipinstalled --ncpus -1 \
    BiocManager \
    shiny \
    shinythemes \
    Seurat\
    markdown\
    shinybusy\
    tidyverse \
    enrichR\
    imager\
    waiter\
    DT\
    igraph\
    readxl\
    shinyWidgets\
    shinyjs\
    RColorBrewer\
    devtools\
    tidyverse\
    magrittr\
    ggrepel\
    data.table\
    matrixStats\
    Matrix

## FROM GITHUB
RUN R -e "devtools::install_github('sqjin/CellChat')"
RUN R -e "devtools::install_github('satijalab/seurat-data')"
RUN R -e "devtools::install_github('mojaveazure/seurat-disk')"
RUN R -e "devtools::install_github('welch-lab/liger')"

## FROM BIOCONDUCTOR
RUN R -e "BiocManager::install('BiocGenerics')"
RUN R -e "BiocManager::install('DelayedArray')"
RUN R -e "BiocManager::install('DelayedMatrixStats')"
RUN R -e "BiocManager::install('limma')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "BiocManager::install('SummarizedExperiment')"
RUN R -e "BiocManager::install('batchelor')"
RUN R -e "BiocManager::install('Matrix.utils')"
RUN R -e "BiocManager::install('rDGIdb')"

RUN R -e "BiocManager::install('SingleR')"
RUN R -e "BiocManager::install('celldex')"
RUN R -e "BiocManager::install('bayNorm')"

#### SPECIFIC R PACKAGES
## MONOCLE3
RUN R -e "devtools::install_github('cole-trapnell-lab/leidenbase')"
RUN R -e "devtools::install_github('cole-trapnell-lab/monocle3')"

## ARCHR
RUN R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
RUN R "ArchR::installExtraPackages()"

## CIPR
RUN R -e "devtools::install_github('atakanekiz/CIPR-Package', build_vignettes = TRUE)"

## SYMPHONY
RUN R -e "devtools::install_github('immunogenomics/symphony')"

## DROPLETUTIL
RUN R -e "BiocManager::install('DropletUtils')"

## ENSEMBLE
RUN R -e "BiocManager::install('ensembldb')"

## SingleCellExperiment (scater)
RUN R -e "BiocManager::install('scater')"


# Copy configuration files into the Docker image
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
RUN rm -rf /srv/shiny-server/*
# COPY /* /srv/shiny-server/

# Get the app code
RUN git clone https://github.com/MillerLab-CPHG/PlaqView.git
RUN cp -r PlaqView/* /srv/shiny-server/

# Make the ShinyApp available at port 80
EXPOSE 80
WORKDIR /srv/shiny-server
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');shiny::runApp('app.R')"

#RUN chown shiny.shiny /usr/bin/shiny-server.sh && chmod 755 /usr/bin/shiny-server.sh

# Run the server setup script
#CMD ["/usr/bin/shiny-server.sh"]
