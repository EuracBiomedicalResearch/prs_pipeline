# Snakemake workflow: prs_pipeline

This pipeline performs the computation of Polygenic Risck Score using different methods.
It integrates different algorithm and allows an easy comparison between the performance of each PRS.


## Installation

An installation of `snakemake` is needed to run this workflow. We suggest to use a `conda`/`mamba` environments to handle the installation, in breif, once `conda` or `mamba` is installed:

```
mamba create -n snakemake -c conda-forge -c bioconda 'snakemake>=8'
```

Activate the conda environment and clone the repository. All the installation and requirements will be handled by `snakemake`.

```
mamba activate snakemake
git clone https://github.com/EuracBiomedicalResearch/prs_pipeline.git
```

## Requirements

In order to compute a Polygenic Risk Score the following files are needed:

1. **Genotype files** in plink (`ped`) format separated by chromosome. For file specification please have a look (here)[https://www.cog-genomics.org/plink/1.9/formats#bed]
2. **GWAS summary statistic**. The results of a GWAS analysis for a specific phenotype
3. Linkage Disequilibrium panel (optional - only for specific methods). The choice of the panel depends on the GWAS analysis. The populations
included in the LD panel should be genetically as close as possible to the population used in the GWAS analysis.

Dataset details can be defined through the configuration file in `YAML` format and `json` files.

Additional information on the configuration can be found in the [`config`](config/README.md).

### Genotype dataset
Genotype dataset can be described through a `json` file, see the example [config/genotype.json](config/genotype.json).

The `json` file should contain the following entries:

- `build`: containing the genome build the samples are genotyped. It could on of `hg19` or `hg38` (`hg18` not supported) as string character.
- `nchrom`: number of chromosomes the genotype files are divided.
- `plinkfiles`: it will be a `python` dictionary with keys starting from 1 to `nchrom`.

### GWAS manifest file

GWAS analysis can be done using a variety of methods, and can be obtained from custom analysis or from 
databases such as [GWAS catalog](https://www.ebi.ac.uk/gwas/).
To support different file specification we use a manifest file in `json` format. Multiple gwases can 
be included in the same manifest file, then the PRS will be computed on each GWAS provided.
You can find and example of this file in the `config` directory:
[config/gwas.json](config/gwas.json) and in the `schema` file
[here](workflow/schemas/gwas_manifest.schema.json)

For each trait, you have to specify the trait type if it is a quantitative or binary trait, and the
columns specification as a python dictionary mapping.

## Available algorithms:

- **PRScs**, see [https://github.com/getian107/PRScs](https://github.com/getian107/PRScs)
- **LDPred2** (through `R` implementation in [`bigsnpr` package](https://privefl.github.io/bigsnpr))
- **lassosum** (implemented in [`bigsnpr` package](https://privefl.github.io/bigsnpr))
- **CT method** (implemented in [`bigsnpr` package](https://privefl.github.io/bigsnpr))
- **sBayesR** (implemented in [gctb](https://cnsgenomics.com/software/gctb/#Overview)) - under
    development

# Running

Once edited the `config.yaml` file you can run the pipeline with the following command.
Needed softwares, packages and files are handled by `snakemake`.

```
cd prs_pipeline
snakemake --cores all --use-conda
```
