#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(rmio)

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
gwasfile <- snakemake@input[["gwas_data"]]

pred_rds <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]

cat(glue("Genotype file: {genofiles}"), "\n")
cat(glue("Mapfile: {mapfile}"), "\n")
cat(glue("gwasfile: {gwasfile}"), "\n")

#---- Setup ----
# Multi-cpu computing
NCORES <- as.integer(snakemake@threads)
cat(glue("Computation on {NCORES} cpus.\n"))

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

#---- Load GWAS ----
gwas_data <- readRDS(gwasfile)
gwas_data <- gwas_data[chr==unique(map_snp$chr)[1]]

#---- Match SNPs ----
# This step has already been done 
# (it's just a double checking that all the SNPs in the bim file are 
# already in the GWAS Data)
# info_snp <- readRDS(snakemake@input[["dfbeta"]])
# TODO: Parametrize `match.min.prop`\
info_snp <- tryCatch(
  snp_match(gwas_data, map_snp, match.min.prop = 0.01),
  error = function(e){ return(NA)}

)

if (is.null(nrow(info_snp))){
  pred_mat <- data.frame(familyID=geno_data$fam$family.ID, sampleID=geno_data$fam$sample.ID, PRS=0)

  saveRDS(pred_mat, file=pred_rds)
  write.table(pred_mat, file=pred_csv, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

} else{
# Check if there is no effects to estimate
# if (length(is_bad) == 0){
is_bad <- rep(FALSE, nrow(info_snp))
# }

#---- Create Beta and LogP vector ----
cat("Creating betas and logs...", "\n")
beta <- rep(NA, ncol(G))
beta[info_snp[!is_bad,]$`_NUM_ID_`] <- info_snp[!is_bad,]$beta
lpval <- rep(NA, ncol(G))
lpval[info_snp[!is_bad,]$`_NUM_ID_`] <- -log10(info_snp[!is_bad,]$p)

t1 <- Sys.time()
#---- Clumping optimizer ----
cat("Start clumping...\n")
clump_res <- snp_grid_clumping(
  G, 
  infos.chr = map_snp$chr,
  infos.pos = map_snp$pos,
  lpS = lpval,
  exclude = which(is.na(lpval)),
  ncores = NCORES)
  #ind.row = ixtrain)
# saveRDS(clump_res, file=clump_res_file)
t2 <- Sys.time()
cat("Clumping optimiziation done in ", t2 - t1, " sec.\n")

t1 <- Sys.time()
#---- Thresholding ----
# TODO: add pvalue threshold as parameter
cat("Start thresholding...\n")
multi_PRS <- snp_grid_PRS(
  G,
  all_keep = clump_res,
  betas = beta,
  lpS = lpval,
  # backingfile = sub(".rds", "", multi_prs_file),
  n_thr_lpS = 10,
  ncores = NCORES)
t2 <- Sys.time()
cat("Thresholding optimiziation done in ", t2 - t1, " sec.\n")

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

# TODO: Check the `nparams`, if less than 1 should report anyway a PRS set to 0
pred_mat <- matrix(0, ncol = nparams, nrow = nrow(multi_PRS))
for (i in seq(1, nparams)){
  ind.col <- seq(i, length.out = 1, by = nparams)
  # pred_mat[, i] <- rowSums(multi_PRS[, ind.col])
  pred_mat[, i] <- multi_PRS[, ind.col]
}
pred_mat <- as.data.frame(pred_mat)
pred_mat <- cbind(familyID=geno_data$fam$family.ID,
sampleID=geno_data$fam$sample.ID, pred_mat)

saveRDS(pred_mat, file=pred_rds)

write.table(pred_mat, file=pred_csv, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
t1 <- Sys.time()
}

#---- Compute PRS using model ----
cat("Quitting...\n")
q(save = "no")

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
