#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(data.table)

#---- Input/Output ----
beta_file <- snakemake@input[["beta_lassosum"]]
df_beta_file <- snakemake@input[["df_beta"]]
geno_file <- snakemake@input[["genotype_rds"]]
pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
map_file <- snakemake@output[["map_file"]]

#---- Read files ----
beta_lasso <- readRDS(beta_file)
geno <- snp_attach(geno_file)
df_beta_lasso <- readRDS(df_beta_file)

#---- Compute prediction ----
pred_lasso <- big_prodMat(geno$genotypes, beta_lasso, ind.col = df_beta_lasso$`_NUM_ID_`)

#---- Save prediction files ----
saveRDS(pred_lasso, file=pred_file)

pred_tab <- data.table(familyID=geno$fam$family.ID, sampleID=geno$fam$sample.ID, prs_lassosum2=pred_lasso)
write.table(pred_tab, file=pred_csv, row.names = FALSE, col.names = TRUE, sep="\t")

#---- Save map file ----
saveRDS(df_beta_lasso, file=map_file)