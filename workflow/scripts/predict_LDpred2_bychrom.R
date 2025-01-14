#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(data.table)

# Avoid nested parallel computation...
Sys.setenv(OPENBLAS_NUM_THREADS=1)

#---- Input/Output ----
# multi_auto_file <- "results/final_mod_ldpred2.rds"
# multi_auto_file <- snakemake@input[[]]
# ytrain_file <- "data/train_pheno_file.rds"
# geno_file <- "data/geno/rds/qc_geno_all.rds"
# df_beta_file <- "data/df_beta_good_with_ld.rds"

# ytrain_file <- snakemake@input[["train_data"]]
# df_beta_file <- snakemake@input[["map_good"]]
# geno_file <- snakemake@input[["genotype_rds"]]
beta_auto_file <- snakemake@input[["pred_auto"]]

pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
map_out <- snakemake@output[["map_file"]]

#---- Resources ----
NCORES <- snakemake@threads

#---- Convert to numeric ----
# This function converts to numeric and set NA to 0
to_double <- function(x){
    y <- as.numeric(x)
    y[!is.na(y)] <- 0
    return(y)
}

#---- Cycle over all files and sum prs up ----
i <- 1
for (rr in beta_auto_file){
    df <- readRDS(rr)
    if (i == 1){
        prs <- to_double(df$PRS)
    } else {
        prs <- prs + to_double(df$PRS)
    }
    i <- i + 1
}

pred_auto <- df[, -PRS]
pred_auto$PRS <- prs

#---- Save predictions into a rdsfile ----
cat("Saving prediction to file rds...\n")
saveRDS(pred_auto, file=pred_file)

#---- Save predictions into a csv file ----
cat("Saving prediction to file csv...\n")
# prs <- data.table(family.ID=fam_info$family.ID, sampleID=fam_info$sample.ID, 
# prs_ldpred2_auto=pred_auto)
write.table(pred_auto, file=pred_csv, sep="\t", 
row.names=FALSE, col.names=FALSE)

#---- Save betas ----
# cat("Saving betas..\n")
# df_beta$beta_prs <- beta_auto
# saveRDS(df_beta, file=map_out)
