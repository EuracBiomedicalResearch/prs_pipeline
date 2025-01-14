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

#---- Input/Ouput ----
# Input
map_rds <- snakemake@input[["map_rds"]]
gwas_rds <- snakemake@input[["gwas_rds"]]
map_ld_rds <- snakemake@input[["map_ld_rds"]]
corfiles <- snakemake@input[["corfiles"]]
geno_file <- snakemake@input[["genotype_rds"]]
pheno <- snakemake@wildcards[["pheno"]]
genotype_conf <- snakemake@params[["genotype_conf"]]
chrom <- snakemake@wildcards[["chrom"]]

# Output
beta_file <- snakemake@output[["beta_file"]]
df_beta_good_rds <- snakemake@output[["df_beta"]]
params_grid_file <- snakemake@output[["params_file"]]
pred_file <- snakemake@output[["pred_file"]]

#---- Resources ----
NCORES <- snakemake@threads
tmpdir <- snakemake@resources[["tmpdir"]]

#---- Load Genotypes ----
geno <- snp_attach(geno_file)
G <- geno$genotypes

#---- Load map_with frequency ----
# map_rds <- "plinkFiles/chrall_map.rds"
map <- readRDS(map_rds)

#---- Load gwas ----
# Remove following lines for local run/debug
# df_beta_rds <- "data/merge_df_beta.rds"
gwas <- readRDS(gwas_rds)

#---- Merge with local genotypes ----
cat("Merging snps...\n")
df_beta <- snp_match(gwas, map)
names(df_beta)[which(names(df_beta) == "_NUM_ID_.ss")] <- "_NUM_ID_.SUMSTAT"

#---- Load training dataset phenotypes ----
# Remove following lines for local run/debug
# traindata_rds <- "data/train_pheno_file.rds"

#---- Load LDref file ----
mapLDref <- readRDS(map_ld_rds)

# Handle hg38 LDreference from bigsnpr
if (genotype_conf$build == "hg38"){
  setnames(mapLDref, "pos", "pos_hg37")
  setnames(mapLDref, "pos_hg38", "pos")
}

#---- Match with LDref ----
df_beta_good <- snp_match(df_beta, mapLDref, match.min.prop = 0)
df_beta_good <- as.data.table(df_beta_good)

#---- Load ld_ref panel ----
## Remove following lines for local run/debug
# ld_ref_rds <- "ld_ref/map_ldref.rds"
# map_ldref <- readRDS(ld_ref_rds)
# Add LD to df_beta
# df_beta[, ld:=map_ldref$ld]

#---- Estimate phenotype variation ----
# TODO: check whether the phenotype is binary or quantitative
# TODO: Update using LD_reference panel frequency
# sd_ldref <- sqrt(2* df_beta$freq * (1 - df_beta$freq))
# sd_ldref <- df_beta[, .(sqrt(2 * freq * (1 - freq)))]

# Quantitative
# TODO: review the current formula (should get an estimation from the original trait used in the GWAS)
# TODO: this is too much dependend by the name of the phenotype, should be standardize
# 
# Estimate trait standard deviation as the first percentile of the GWAS effects
# See: https://privefl.github.io/bigsnpr/articles/LDpred2.html for further details
# sd_trait_est <- quantile(df_beta[, .(sqrt(0.5 * (n_eff * beta_se**2) + beta**2))]$V1, probs=0.01)

# sd_train_est <- sqrt(0.016)
# sd_ss <- df_beta[, .(1 / sqrt(n_eff * beta_se**2 + beta**2))]
# sd_ss <- as.numeric(sd_ss$V1)
# sd_ss <- sd_ss / quantile(sd_ss, 0.99) * sqrt(0.5)

# Binary case
# sd_ss <- df_beta[, .(2 / sqrt(n_eff * beta_se**2 + beta**2))]

#---- Mark bad snps ----
#sd_ss < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
# sd_ss < 0.1 | sd_af < 0.05)
#
# is_bad <- 
#   sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05 #| df_beta$freq < 0.01 | df_beta$freq > 0.99

#---- Subselect only GOOD snps ----
# cat("Select only good snps...\n")
# df_beta_good <- df_beta[!is_bad, ]

#---- Remove unused files ----
# cat("Removing unused files...and garbage collecting...\n")
# rm(df_beta, map)
# gc()

#---- Create LD column on df_beta ----
# df_beta_good[, ld:=0]

#---- Read correlation matrix previously computed ----
# Remove the following line for local execution
# tmpdir <- "tmp-data"
tmp <- tempfile(tmpdir = tmpdir)
corr <- NULL
t1 <- Sys.time()
cat(chrom, ".. ", sep = "")
ix <- df_beta_good[chr == chrom, `_NUM_ID_`]

# retrieve indexes relative to the correlation matrix positions
ix2 <- match(ix, which(mapLDref$chr == chrom))

# Get correlation file
# fi <- grepl(glue("chr{mychr}.rds"), corfiles)
ff <- corfiles

cat("read matrix...\n")
corr_chr <- readRDS(ff)[ix2, ix2]

corr <- as_SFBM(corr_chr, tmp, compact = TRUE)

t2 <- Sys.time()
cat(glue("Compute correlation in {format(t2-t1)}...\n"))

#---- Compute lassosum2 betas ----
cat("Computing beta with lassosum2...")
t1 <- Sys.time()
beta_lassosum2 <- snp_lassosum2(corr, df_beta_good, ncores = NCORES)
t2 <- Sys.time()
cat(glue("in {format(t2 - t1)}"), "\n Saving...\n")
saveRDS(beta_lassosum2, file=beta_file)

#---- Lassosum2 parameters ----
params2 <- attr(beta_lassosum2, "grid_param")
saveRDS(params2, file=params_grid_file)

#---- Compute predictions ----
setnames(df_beta_good, "_NUM_ID_.ss", "_NUM_ID_.geno_chrom")
pred <- big_prodMat(G, beta_lassosum2, ind.col = df_beta_good[["_NUM_ID_.geno_chrom"]])
colnames(pred) <- paste0("p", seq(1, nrow(params2)))

#---- Write prediction ----
cat("Saving prediction for lassosum2\n")
pred_grid_df <- data.table(geno$fam)
pred_grid_df <- cbind(pred_grid_df, pred)
saveRDS(pred_grid_df, file=pred_file)

#---- Export df_beta_good for future computation ----
cat("Saving df_beta_good\n")
saveRDS(df_beta_good, file=df_beta_good_rds)