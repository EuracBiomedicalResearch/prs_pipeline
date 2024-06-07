# Configuration

Different option can be change in the [config/config.yaml](config/config.yaml) file.

## Output

With this set of options it's possible to change the output results folder and download resources
path. All paths are **relative** to the running location. Either absolute path are accepted.

Results will be in the directory specified by `output_dir`, while quality control on genotype will
be in the `data_dir`.

```
output_dir: "results"
data_dir: "data"
```

If some resources are already available from previous computation it's possible to reuse them
specifing the `cache_dir` parameter.

```
# Cache directory for resources:
cache_dir: "/storage03/mfilosi/CKDgenPRS/prs_pipeline/resources"
```

## Genotype QC

Parameters to use for quality controll process of the genotype. These are passed directly to `plink`
for quality control with the relative options.
See [plink](https://www.cog-genomics.org/plink/1.9/index) for the meaning of the options.

```
# Genotype QC parameters:
maf: 0.05
hwe: 1e-50
genorate: 0.1
mind: 0.1
```

## LD reference

Definition of the reference panel for PRScs algorithm.
Available options for population and data for LD panels are available
[here](https://github.com/getian107/PRScs). Needed file are downloaded automatically based on the 
parameters defined here:

```
ld_reference:
  population: "eur"
  data: "1KG"
```

# Algorithms

Activation of the algorithm to use for the analysis.

```
# PRS-method activation
prs_algorithms:
  ldpred2:
    activate: true
  lassosum2:
    activate: true
  prscs:
    activate: true
  sct:
    activate: true
```
