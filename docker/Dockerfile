# R environment for the RNA fingerprinting pipeline.
# Build LOCALLY, push to Docker Hub, then pull/convert to Apptainer on Trillium
# (compute nodes have no internet; all installs happen here at build time).
#
#   docker build -t <dockerhub-user>/rna-fingerprinting:latest -f docker/Dockerfile .
#   docker push  <dockerhub-user>/rna-fingerprinting:latest
#   # on Trillium login node:
#   apptainer pull rna-fingerprinting.sif docker://<dockerhub-user>/rna-fingerprinting:latest
#
# Seurat is preinstalled in the base image; we add the Satija RNA fingerprinting
# package + its precomputed-dictionary data package.

FROM satijalab/seurat:5.0.0

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('remotes','stringr'), repos='https://cloud.r-project.org')"

# RNA fingerprinting method + precomputed dictionaries (GWPS, Immune Dictionary).
RUN R -e "remotes::install_github('satijalab/rna-fingerprinting', upgrade='never')"
RUN R -e "remotes::install_url('https://github.com/satijalab/rna-fingerprinting-data/raw/refs/heads/main/RNAFingerprintingData_0.1.0.tar.gz')"

# Sanity check that the package loads.
RUN R -e "library(RNAFingerprinting); cat('RNAFingerprinting OK\n')"

WORKDIR /work
