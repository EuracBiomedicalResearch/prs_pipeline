FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="b03305cb9fee89d850d0ff8554fe29c112dff5bda687876602e3196e2fba66a0"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/bigsnpr.yaml
#   prefix: /conda-envs/5c0e31c1f241bef4391aa8e9b73d1f11
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
RUN mkdir -p /conda-envs/5c0e31c1f241bef4391aa8e9b73d1f11
COPY workflow/envs/bigsnpr.yaml /conda-envs/5c0e31c1f241bef4391aa8e9b73d1f11/environment.yaml

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
#   prefix: /conda-envs/b5781f6e09c9a885bd5f0205307ed8b4
#   channels:
#     - defaults
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - prscs=1.1.0
#     - scipy=1.11
#     - h5py
#     - numpy
RUN mkdir -p /conda-envs/b5781f6e09c9a885bd5f0205307ed8b4
COPY workflow/envs/prscs.yaml /conda-envs/b5781f6e09c9a885bd5f0205307ed8b4/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5c0e31c1f241bef4391aa8e9b73d1f11 --file /conda-envs/5c0e31c1f241bef4391aa8e9b73d1f11/environment.yaml && \
    mamba env create --prefix /conda-envs/fab2051204386f30216706b050891405 --file /conda-envs/fab2051204386f30216706b050891405/environment.yaml && \
    mamba env create --prefix /conda-envs/b5781f6e09c9a885bd5f0205307ed8b4 --file /conda-envs/b5781f6e09c9a885bd5f0205307ed8b4/environment.yaml && \
    mamba clean --all -y
