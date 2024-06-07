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
df_beta_file <- snakemake@input[["map_good"]]
geno_file <- snakemake@input[["genotype_rds"]]
beta_auto_file <- snakemake@input[["beta_auto"]]

pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
map_out <- snakemake@output[["map_file"]]

#---- Resources ----
NCORES <- snakemake@threads

#---- Read RDS from the input ----
# multi_auto <- readRDS(multi_auto_file)
# ytrain <- readRDS(ytrain_file)
df_beta <- readRDS(df_beta_file)
geno <- snp_attach(geno_file)
beta_auto <- readRDS(beta_auto_file)

#---- Fam info from genotype ----
fam_info <- geno$fam
# fam_info <- fam_info %>% mutate(
#   sample.ID=str_pad(sample.ID, width=10, side="left", pad=0)
# )

#---- Subset train data and genotype data ----
# ixsamp <- match(fam_info$sample.ID, ytrain$AID)

# NB create `ixtrain` object that contains the samples ID corresponding the the
# genotype file for the sample ID available in the phenotype
# If all samples are found in the phenotype and genotype `ixtrain` is equal to 
# all rows in the genotype file, otherwise will be subsetted.
# The model fitting require the phenotypes to be **already** subsetted matching 
# the genotype sample order!
# if (sum(is.na(ixsamp)) > 0){
#   ixtrain <- which(!is.na(ixsamp))
#   ytrain_sub <- ytrain[ixsamp[!is.na(ixsamp)],]
# } else {
#   ixtrain <- 1:length(ixsamp)
#   ytrain_sub <- ytrain[ixsamp,]
# }

#---- Compute PRS from LDpred2 ----
G <- geno$genotypes
cat("Computing predictions...\n")
pred_auto <- big_prodVec(G, beta_auto, 
                         ind.col = df_beta$`_NUM_ID_`, 
                         ncores=NCORES)

#---- Save predictions into a rdsfile ----
cat("Saving prediction to file rds...\n")
saveRDS(pred_auto, file=pred_file)

#---- Save predictions into a csv file ----
cat("Saving prediction to file csv...\n")
prs <- data.table(family.ID=faminfo$family.ID, sampleID=faminfo$sample.ID, prs_ldpred2_auto=pred_auto)
write.table(prs, file=pred_csv, sep="\t", row.names=FALSE, col.names=FALSE)

#---- Save betas ----
cat("Saving betas..\n")
df_beta$beta_prs <- beta_auto
saveRDS(df_beta, file=map_out)
