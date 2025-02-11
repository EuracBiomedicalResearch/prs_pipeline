FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="3baaf85514028be1d7a542a762f2a5eee8cc3a5f1a19ceaeaa73c22cc9093714"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bcftools.yaml
#   prefix: /conda-envs/2c0d46c7ebcc3984dfdd9668f1a3abba
#   channels:
#     - defaults
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools
#     - bcftools
#     - tabix
RUN mkdir -p /conda-envs/2c0d46c7ebcc3984dfdd9668f1a3abba
COPY workflow/envs/bcftools.yaml /conda-envs/2c0d46c7ebcc3984dfdd9668f1a3abba/environment.yaml

# Conda environment:
#   source: workflow/envs/bigsnpr.yaml
#   prefix: /conda-envs/124bf7e9aaf9d2192c1ea1b7d35c4a6b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - r-stringr =1.5
#     - r-bigsnpr =1.12.2
#     - r-bigreadr =0.2.5
#     - r-forcats
#     - r-data.table
#     - r-ggplot2
#     - r-glue
#     - r-rmio
#     - r-jsonlite
#     - r-dplyr
#     - r-R.utils
#     - r-tidyr
RUN mkdir -p /conda-envs/124bf7e9aaf9d2192c1ea1b7d35c4a6b
COPY workflow/envs/bigsnpr.yaml /conda-envs/124bf7e9aaf9d2192c1ea1b7d35c4a6b/environment.yaml

# Conda environment:
#   source: workflow/envs/plink.yaml
#   prefix: /conda-envs/fab2051204386f30216706b050891405
#   channels:
#     - defaults
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - plink =1.90b6.21
RUN mkdir -p /conda-envs/fab2051204386f30216706b050891405
COPY workflow/envs/plink.yaml /conda-envs/fab2051204386f30216706b050891405/environment.yaml

# Conda environment:
#   source: workflow/envs/prscs.yaml
#   prefix: /conda-envs/6378c85b56c9bb25d29b90c8204b496c
#   channels:
#     - defaults
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     # - prscs=1.1.0
#     - scipy=1.11
#     - h5py
#     - numpy
#     - liftover
#     - pip
#     - pip:
#       - git+https://github.com/filosi/PRScs.git@feature/package
RUN mkdir -p /conda-envs/6378c85b56c9bb25d29b90c8204b496c
COPY workflow/envs/prscs.yaml /conda-envs/6378c85b56c9bb25d29b90c8204b496c/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/2c0d46c7ebcc3984dfdd9668f1a3abba --file /conda-envs/2c0d46c7ebcc3984dfdd9668f1a3abba/environment.yaml && \
    mamba env create --prefix /conda-envs/124bf7e9aaf9d2192c1ea1b7d35c4a6b --file /conda-envs/124bf7e9aaf9d2192c1ea1b7d35c4a6b/environment.yaml && \
    mamba env create --prefix /conda-envs/fab2051204386f30216706b050891405 --file /conda-envs/fab2051204386f30216706b050891405/environment.yaml && \
    mamba env create --prefix /conda-envs/6378c85b56c9bb25d29b90c8204b496c --file /conda-envs/6378c85b56c9bb25d29b90c8204b496c/environment.yaml && \
    mamba clean --all -y
