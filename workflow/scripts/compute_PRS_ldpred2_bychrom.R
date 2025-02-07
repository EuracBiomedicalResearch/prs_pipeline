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
# df_beta_rds <- snakemake@input[["df_beta"]]
# traindata_rds <- snakemake@input[["train_data"]]
map_rds <- snakemake@input[["map_rds"]]
gwas_rds <- snakemake@input[["gwas_rds"]]
map_ld_rds <- snakemake@input[["map_ld_rds"]]
# ld_ref_rds <- snakemake@input[["ld_ref"]]
heritability_file <- snakemake@input[["heritability"]]
corfiles <- snakemake@input[["corfiles"]]
geno_file <- snakemake@input[["genotype_rds"]]
pheno <- snakemake@wildcards[["pheno"]]
genotype_conf <- snakemake@params[["genotype_conf"]]
chrom <- snakemake@wildcards[["chrom"]]

# Output
model_file <- snakemake@output[["final_mod"]]
beta_file <- snakemake@output[["beta_auto"]]
df_beta_good_rds <- snakemake@output[["df_good"]]
beta_grid_file <- snakemake@output[["beta_grid"]]
pred_grid_file <- snakemake@output[["pred_grid"]]
params_grid_file <- snakemake@output[["params_grid"]]
pred_auto_file <- snakemake@output[["pred_auto"]]

#---- Resources ----
NCORES <- snakemake@threads
tmpdir <- snakemake@resources[["tmpdir"]]

#---- Load Genotypes ----
geno <- snp_attach(geno_file)
G <- geno$genotypes

#---- Load map_with frequency ----
# map_rds <- "plinkFiles/chrall_map.rds"
map <- readRDS(map_rds)
map <- geno$map
names(map)[-3] <- c("chr", "rsid", "pos", "a1", "a0")

#---- Load gwas ----
gwas <- readRDS(gwas_rds)
gwas <- gwas[chr==chrom,]

#---- Merge with local genotypes ----
cat("Merging snps...\n")
df_beta <- tryCatch(
  snp_match(gwas, map, match.min.prop = 0.01),
  error = function(e){ return(NA)}
)

if (is.null(nrow(df_beta))){
  pred_mat <- data.frame(familyID=geno$fam$family.ID, sampleID=geno$fam$sample.ID, PRS=0)
  saveRDS(pred_mat, file=pred_auto_file)
  saveRDS(pred_mat, file=pred_grid_file)

  # Touch files
  file.create(df_beta_good_rds)
  file.create(params_grid_file)
} else {
  # Rename df_beta columns
  names(df_beta)[which(names(df_beta) == "_NUM_ID_.ss")] <- "_NUM_ID_.SUMSTAT"

  #---- Load LDref file ----
  cat(glue("Load LD ref file {map_ld_rds}"), "\n")
  t1 <- Sys.time()
  mapLDref <- readRDS(map_ld_rds)
  cat(glue("Loaded in {Sys.time() - t1}"), "\n")

  # Handle hg38 LDreference from bigsnpr
  if (genotype_conf$build == "hg38"){
    setnames(mapLDref, "pos", "pos_hg37")
    setnames(mapLDref, "pos_hg38", "pos")
  }

  #---- Match with LDref ----
  # TODO: add check on this step if there are no matching variants
  df_beta_good <- snp_match(df_beta, mapLDref, match.min.prop = 0)
  df_beta_good <- as.data.table(df_beta_good)

  #---- Read correlation matrix previously computed ----
  tmp <- tempfile(tmpdir = tmpdir)
  ix <- df_beta_good[chr == chrom, `_NUM_ID_`]

  # retrieve indexes relative to the correlation matrix positions
  ix2 <- match(ix, which(mapLDref$chr == chrom))

  # Get correlation file
  cat(glue("Reading {corfiles}"), "\n")
  corr_chr <- readRDS(corfiles)[ix2, ix2]
  corr <- as_SFBM(corr_chr, tmp, compact = TRUE)

  #---- Export df_beta_good for future computation ----
  saveRDS(df_beta_good, file=df_beta_good_rds)

  #---- Heritability estimation of LD score regression ----
  # to be used as a starting value in LDpred2-auto
  ldsc <- readRDS(heritability_file)
  h2_est <- ldsc[["h2"]]

  # LDpred2-auto
  t1 <- Sys.time()
  cat("Running LDpred-auto estimation...\n")
  multi_auto <- snp_ldpred2_auto(corr, df_beta_good, h2_init = h2_est,
                                vec_p_init = seq_log(1e-4, 0.9, 
                                length.out = 10),
                                allow_jump_sign = FALSE, shrink_corr = 0.95,
                                ncores = NCORES)
  t2 <- Sys.time()
  cat(glue("LDpred model run in {format(t2-t1)}"),"\n")
  cat(glue("Saving {model_file}..."), "\n")

  #---- Extract parameters from multi-auto model ----
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

  if (length(keep) == 0){
    keep <- 1:length(multi_auto)
  }
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  df_beta_good[, beta_auto:=beta_auto]
  setnames(df_beta_good, "_NUM_ID_.ss", "_NUM_ID_.geno_chrom")

  #---- Predict PRS using beta auto model -----
  pred_auto <- big_prodVec(G, beta_auto,
    ind.col = df_beta_good[["_NUM_ID_.geno_chrom"]])
  pred_df <- data.table(geno$fam)
  pred_df[, PRS := pred_auto]

  # Saving...
  saveRDS(pred_df, file=pred_auto_file)

  #---- LDpred2 - grid ----
  # Running grid model for LDPRED2.
  cat("Running ldpred2-grid model\n")
  # Define grid parameters
  # TODO: add this as customizable in the settings
  h2_seq <- round(h2_est * c(0.01, 0.3, 0.7, 1, 1.4, 2), 4)
  h2_seq <- h2_seq[h2_seq > 0]
  p_seq <- signif(seq_log(1e-5, 1, length.out = 10), 2)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
  # Saving grid parameters
  saveRDS(params, file=params_grid_file)

  # Compute the grid model weight
  beta_grid <- snp_ldpred2_grid(corr, df_beta_good, params, ncores = NCORES)
  # saveRDS(beta_grid, file=beta_grid_file)

  # Compute the prediction
  pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta_good[["_NUM_ID_.geno_chrom"]])
  colnames(pred_grid) <- paste0("p", seq(1, nrow(params)))
  pred_grid_df <- data.table(geno$fam)
  pred_grid_df <- cbind(pred_grid_df, pred_grid)
  saveRDS(pred_grid_df, file=pred_grid_file)
}

#--- Quit ----
quit(save = "no")