#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)#, lib.loc = "/home/mfilosi/R/rocker-rstudio/4.2")
library(rmio)

#/shared/bioinf/R/bin/Rscript-4.3-BioC3.17

# Avoid nested parallel computation...
Sys.setenv(OPENBLAS_NUM_THREADS=1)

#---- Parameter ----
#| force: if TRUE then ignore existing files and recompute the model from scratch
force <- snakemake@params[["force"]]
force <- as.logical(force)

cat(glue('Force parameter {force}'), "\n")

#---- Input/Output ----
genofiles <- snakemake@input[["geno_data"]]
mapfile <- snakemake@input[["map_data"]]
trainfile <- snakemake@input[["train_data"]]
gwasfile <- snakemake@input[["gwas_data"]]

pred_rds <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]

cat(glue("Genotype file: {genofiles}"), "\n")
cat(glue("Mapfile: {mapfile}"), "\n")
cat(glue("trainfile: {trainfile}"), "\n")
cat(glue("gwasfile: {gwasfile}"), "\n")


#---- Setup ----
# proj_path <- "/scratch/mfilosi/CKDgenPRS"
# data_path <- file.path(proj_path, "data")
# geno_path <- file.path(proj_path, "plinkFiles")
# tmp_path <- file.path(proj_path, "tmp-data")

# Multi-cpu computing
NCORES <- as.integer(snakemake@threads)
cat(glue("Computation on {NCORES} cpus.\n"))

#---- Get phenotypes and training data ----
# train_data <- readRDS(file.path(data_path, "train_pheno_file.rds"))
# train_data <- readRDS(trainfile)

t1 <- Sys.time()
#---- Load Genotypes ----
cat("Reading genotype data...\n")
cat(glue("genotype file: {genofiles}"), "\n")
geno_data <- snp_attach(genofiles)
# Get time and print log
t2 <- Sys.time()
cat("Genotype data readed in ", t2 - t1, " sec.\n")

#---- Genotype info ----
G <- geno_data$genotypes

#---- Map info ----
map_snp <- readRDS(mapfile)
# names(map_snp) <- c("chr", "rsid", "genetic_dist", "pos", "a0", "a1")

cat(glue("Read {ncol(G)} snps from genotype and {nrow(map_snp)} from {mapfile}") , "\n")

#---- Fam info from genotype ----
fam_info <- geno_data$fam
fam_info <- fam_info %>% mutate(
  sample.ID=str_pad(sample.ID, width=10, side="left", pad=0)
)

#---- Subset train data and genotype data ----
# ixsamp <- match(fam_info$sample.ID, train_data$AID)

# NB create `ixtrain` object that contains the samples ID corresponding the the
# genotype file for the sample ID available in the phenotype
# If all samples are found in the phenotype and genotype `ixtrain` is equal to 
# all rows in the genotype file, otherwise will be subsetted.
# The model fitting require the phenotypes to be **already** subsetted matching 
# the genotype sample order!
# if (sum(is.na(ixsamp)) > 0){
#   ixtrain <- which(!is.na(ixsamp))
#   train_data <- train_data[ixsamp[!is.na(ixsamp)],]
# } else {
#   ixtrain <- 1:length(ixsamp)
#   train_data <- train_data[ixsamp,]
# }

#---- Phenotype to train ----
# y <- train_data %>% dplyr::pull(eGFR1)
# Rescale phenotype
# yscaled <- scale(y, center=TRUE, scale=1)

#---- Create covariates ----
# covariates <- train_data %>% 
#   dplyr::select(sex, age, starts_with("PC"))
# cov_mat <- covar_from_df(covariates)

#---- Load GWAS ----
# gwas_data <- readRDS(file.path(data_path, "19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.rds"))
gwas_data <- readRDS(gwasfile)

#---- Match SNPs ----
# This step has already been done 
# (it's just a double checking that all the SNPs in the bim file are 
# already in the GWAS Data)
# info_snp <- readRDS(snakemake@input[["dfbeta"]])
info_snp <- snp_match(gwas_data, map_snp)

# TODO: Need to add choice if it's binary or quantitative trait
sd_ldref <- sqrt(2* map_snp[info_snp$`_NUM_ID_`, "freq"] * (1 - map_snp[info_snp$`_NUM_ID_`, "freq"]))
sd_trait_est <- quantile(sqrt(0.5 * (info_snp$n_eff * info_snp$beta_se**2) + info_snp$beta**2), probs=0.99)
sd_ss <- 1 / sqrt(info_snp$n_eff * info_snp$beta_se**2 + info_snp$beta**2)
sd_ss <- sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)

is_bad <- 
  sd_ss < (0.7 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

# Check if there is no effects to estimate
if (length(is_bad) == 0){
  is_bad <- rep(FALSE, nrow(info_snp))
}

#---- Create Beta and LogP vector ----
cat("Creating betas and logs...", "\n")
beta <- rep(NA, ncol(G))
beta[info_snp[!is_bad,]$`_NUM_ID_`] <- info_snp[!is_bad,]$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp[!is_bad,]$`_NUM_ID_`] <- -log10(info_snp[!is_bad,]$p_value)

t1 <- Sys.time()
#---- Clumping optimizer ----
clump_res_file <- snakemake@output[["clump_opt"]]
if (file.exists(clump_res_file) & force==FALSE){
  cat("Found clumping file, reading...\n")
  clump_res <- readRDS(clump_res_file)
} else {
  cat("Start clumping...\n")
  clump_res <- snp_grid_clumping(
    G, 
    infos.chr = map_snp$chr, 
    infos.pos = map_snp$pos,
    lpS = lpval,
    exclude = which(is.na(lpval)),
    ncores=NCORES)
    #ind.row = ixtrain)
  saveRDS(clump_res, file=clump_res_file)
  t2 <- Sys.time()
  cat("Clumping optimiziation done in ", t2 - t1, " sec.\n")
}

t1 <- Sys.time()
#---- Thresholding ----
multi_prs_file <- snakemake@output[["multi_prs"]]
multi_prsbk_file <- snakemake@output[["multi_prs_bk"]]
if (file.exists(multi_prs_file) & force == FALSE){
  cat("Found thresholding file, reading...\n")  
  multi_PRS <- readRDS(multi_prs_file)  
} else {
  if (file.exists(multi_prsbk_file)){
    file.remove(multi_prsbk_file)
  }
  cat("Start thresholding...\n")
  multi_PRS <- snp_grid_PRS(
    G, 
    all_keep = clump_res, 
    betas=beta, 
    lpS=lpval, 
    # ind.row = ixtrain,
    backingfile = sub(".rds", "", multi_prs_file), 
    n_thr_lpS = 10, 
    ncores = NCORES)
  t2 <- Sys.time()
  cat("Thresholding optimiziation done in ", t2 - t1, " sec.\n")
}

# Save values into a matrix
# 
# NB multi_PRS contains the PRS for each of the parameters in the grid.ldS.thr and all_keep grid.
# store by chromosome, lpthreshold and grid parameters.
# indices from 1 to nrow(attr(all_keep, "grid")) * length(lpS_thr) are for chromosome 1 and so on
all_keep <- attr(multi_PRS, "all_keep")
lpS_thr <- attr(multi_PRS, "grid.lpS.thr")

grids <- attr(all_keep, "grid")
ngrids <- nrow(grids)
n_thr <- length(lpS_thr)
nparams <- ngrids * n_thr
nchroms <- 22

pred_mat <- matrix(0, ncol=nparams, nrow=nrow(multi_PRS))
for (i in seq(1, nparams)){
  ind.col <- seq(i, length.out=nchroms, by=nparams)
  pred_mat[, i] <- rowSums(multi_PRS[, ind.col])
}
pred_mat <- as.data.frame(pred_mat)
pred_mat <- cbind(familyID=geno_data$fam$family.ID, sampleID=geno_data$fam$sample.ID, pred_mat)
saveRDS(pred_mat, file=pred_rds)
write.table(pred_mat, file=pred_csv, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
t1 <- Sys.time()


#---- Compute PRS using model ----
# final_mod_file <- snakemake@output[["final_mod"]]
# if (file.exists(final_mod_file) & force == FALSE ){
#   final_mod <- readRDS(final_mod_file)
#   cat("Load `final_models` scores....")
#   # cat("Nothing to be done for this script, `final_mod` has been already computed!")
# } else {
#   cat("Start PRS model optimization...\n")
#   final_mod <- snp_grid_stacking(
#     multi_PRS, 
#     yscaled, 
#     # ind.train = ixtrain, 
#     covar.train=cov_mat, 
#     ncores = NCORES, K = 10)
#   t2 <- Sys.time()
#   # Save final model
#   saveRDS(final_mod, file=final_mod_file)
#   cat("Model optimiziation done in ", t2 - t1, " sec.\n")
# }

cat("Quitting...\n")
q(save="no")

#---- Prediction ----
newbeta <- final_mod$beta.G
indbeta <- which(newbeta != 0)
# Add prediction based on covariates
covpred <- cov_mat[ixtrain, ] %*% final_mod$beta.covar
saveRDS(covpred, snakemake@output[["cov_pred"]])
# saveRDS(covpred, file.path(data_path, "covariate_predictions.rds"))

# Add prediction based on genetic only
genpred <- big_prodVec(G, newbeta[indbeta], ind.row = ixtrain, ind.col = indbeta)
saveRDS(genpred, snakemake@output[["gen_pred"]])
# saveRDS(genpred, file.path(data_path, "genetic_predictions.rds"))

#---- Add all prediction ----
pred <- final_mod$intercept + genpred + covpred
saveRDS(covpred, snakemake@output[["pred"]])
# saveRDS(pred, file.path(data_path, "predictions.rds"))

#---- Quitting ----
cat("Quitting...\n")
q(save="no")

#---- Regression measures ----
# Root mean squared error
rmse <- function(y, ypred){
  rmse <- sqrt(sum((y - ypred)**2)/ length(y))
  return(rmse)
}

# Coefficient of determination
R2 <- function(y, ypred){
  m <- mean(y)
  sstot <- sum((y - m)**2)
  ssres <- sum((y - ypred)**2)
  r2 <- ssres / sstot
  return(r2)
}

mymod <- bigstatsr::big_spLinReg(G, yscaled, ind.train = ixtrain, covar.train = cov_mat, )
