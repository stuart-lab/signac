FROM satijalab/seurat:latest

# Install system dependencies
RUN apt-get update --fix-missing \
  && apt-get install -y software-properties-common \
  && add-apt-repository -y ppa:git-core/ppa \
  && apt-get update \
  && apt-get install -y git libbz2-dev liblzma-dev

# Install Bioconductor dependencies

RUN R --slave --no-restore --no-save -e "install.packages('BiocManager'); BiocManager::install(c('GenomeInfoDbData', 'GO.db'))"

# Install Signac

RUN  R --slave --no-restore --no-save -e "setRepositories(ind=1:2); install.packages('Signac',  repos = c('https://cloud.r-project.org', 'https://bioconductor.org/packages/3.13/bioc'),dependencies = TRUE)"

CMD [ "R" ]
