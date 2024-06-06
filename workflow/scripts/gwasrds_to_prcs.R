#!/shared/bioinf/R/bin/Rscript-4.3-BioC3.17
#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(data.table)

#---- Input/Output ----
gwas_rds <- snakemake@input[["gwas_rds"]]

#---- Read in the GWAS ----
cat("Reading file...\n")
gwas <- readRDS(gwas_rds)

#---- Extract only needed columns ----
gwas <- gwas[freq > 0.01 & freq < 0.99, c("id", "a0", "a1", "beta", "beta_se")]
names(gwas) <- c("SNP", "A1", "A2", "BETA", "SE")

#---- Write gwas output ----
cat("Writing file...\n")
data.table::fwrite(gwas, file=snakemake@output[["gwas_prscs"]], sep="\t")
