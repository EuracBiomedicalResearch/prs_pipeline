#!/shared/bioinf/R/bin/Rscript-4.3-BioC3.17

#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr, lib.loc = "/home/mfilosi/R/rocker-rstudio/4.2")
library(data.table)

#---- Input/Output ----
beta_file <- snakemake@input[["beta_lassosum"]]
df_beta_file <- snakemake@input[["df_beta"]]
geno_file <- snakemake@input[["genotype_rds"]]
pred_file <- snakemake@output[["pred"]]

#---- Read files ----
beta_lasso <- readRDS(beta_file)
geno <- snp_attach(geno_file)
df_beta_lasso <- readRDS(df_beta_file)


#---- Compute prediction ----
pred_lasso <- big_prodMat(geno$genotypes, beta_lasso, ind.col = df_beta_lasso$`_NUM_ID_`)

#---- Save file ----
saveRDS(pred_lasso, file=pred_file)