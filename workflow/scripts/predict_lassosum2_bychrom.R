#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(data.table)

#---- Input/Output ----
# beta_file <- snakemake@input[["beta_lassosum"]]
# df_beta_file <- snakemake@input[["df_beta"]]
# geno_file <- snakemake@input[["genotype_rds"]]
pred_lasso <- snakemake@input[["pred_lassosum"]]
pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
# map_file <- snakemake@output[["map_file"]]

#---- Convert to numeric ----
# This function converts to numeric and set NA to 0
to_double <- function(x){
    y <- as.numeric(x)
    y[is.na(y)] <- 0
    return(y)
}

#---- Cycle over all files and sum prs up ----
i <- 1
for (rr in pred_lasso){
    df <- readRDS(rr)
    if (i == 1){
        prs <- df[, 7:ncol(df)]
    } else {
        prs <- prs + df[, 7:ncol(df)]
    }
    i <- i + 1
}
# Add family ID and sample ID
prs_all_samps <- cbind(df[, 1:7], prs)

#---- Save predictions into a rdsfile ----
cat("Saving prediction to file rds...\n")
saveRDS(prs_all_samps, file = pred_file)

#---- Save predictions into a csv file ----
cat("Saving prediction to file csv...\n")
# prs <- data.table(family.ID=fam_info$family.ID, sampleID=fam_info$sample.ID, 
# prs_ldpred2_auto=pred_auto)
write.table(prs_all_samps, file = pred_csv, sep = "\t", row.names = FALSE,
    col.names = FALSE)
