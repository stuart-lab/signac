FROM satijalab/seurat:latest

# Install system dependencies
RUN apt-get update --fix-missing
RUN apt-get install -y software-properties-common
RUN add-apt-repository -y ppa:git-core/ppa
RUN apt-get update
RUN apt-get install -y git

# Install Bioconductor dependencies

RUN R --slave --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --slave --no-restore --no-save -e "BiocManager::install(c('GenomeInfoDbData', 'GO.db'))"

# Install Signac

RUN  R --slave --no-restore --no-save -e "setRepositories(ind=1:2); install.packages('Signac',  repos = c('https://cloud.r-project.org', 'https://bioconductor.org/packages/3.11/bioc'),dependencies = TRUE)"

CMD [ "R" ]
