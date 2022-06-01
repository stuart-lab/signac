FROM satijalab/seurat:latest

RUN apt-get update --fix-missing \
  && apt-get install -y libbz2-dev liblzma-dev \
  &&  R --slave --no-restore --no-save -e "options(repos=c(CRAN='https://cloud.r-project.org')); setRepositories(ind=1:3); install.packages('Signac', dependencies = TRUE)"

CMD [ "R" ]
