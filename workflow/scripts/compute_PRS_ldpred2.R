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
# TODO: define in/out parameters from snakemake
# df_beta_rds <- snakemake@input[["df_beta"]]
# traindata_rds <- snakemake@input[["train_data"]]
map_rds <- snakemake@input[["map_rds"]]
gwas_rds <- snakemake@input[["gwas_rds"]]
map_ld_rds <- snakemake@input[["map_ld_rds"]]
# ld_ref_rds <- snakemake@input[["ld_ref"]]
ldfiles <- snakemake@input[["ldfiles"]]
corfiles <- snakemake@input[["corfiles"]]
geno_file <- snakemake@input[["genotype_rds"]]
pheno <- snakemake@wildcards[["pheno"]]
genotype_conf <- snakemake@params[["genotype_conf"]]

# Remove following lines for local run/debug
# ldfiles <- list.files("ld_ref", pattern="ld_chr")
# ldfiles <- file.path("ld_ref", ldfiles)
print(ldfiles)

# Output
heritability_file <- snakemake@output[[1]]
model_file <- snakemake@output[[2]]
beta_file <- snakemake@output[[3]]
df_beta_good_rds <- snakemake@output[[4]]
beta_grid_file <- snakemake@output[[5]]
pred_grid_file <- snakemake@output[[6]]
params_grid_file <- snakemake@output[[7]]

print(genotype_conf$build)

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
# df_beta <- as.data.table(df_beta

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

#---- Load training dataset phenotypes ----
# Remove following lines for local run/debug
# traindata_rds <- "data/train_pheno_file.rds"
# traindata <- readRDS(traindata_rds)

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

#---- Plot snp std deviation estimation ----
# TODO: Add plotting if requested from the snakemake pipeline
# df2plot <- data.frame(sd_ldref=sd_ldref$V1, sd_ss=sd_ss$V1,  is_bad=is_bad)
# 
# ggplot(df2plot, aes(sd_ldref, sd_ss, color=is_bad)) + 
#   geom_point(alpha=0.5) + theme_bigstatsr() + coord_equal() + 
#   geom_abline(linetype = 2, color = "red") +
#   labs(x = "Standard deviations derived from allele frequencies of the LD reference",
#        y = "Standard deviations derived from the summary statistics",
#        color = "Removed?")

#---- Subselect only GOOD snps ----
# cat("Select only good snps...\n")
# df_beta_good <- df_beta[!is_bad, ]

#---- Remove unused files ----
# cat("Removing unused files...and garbage collecting...\n")
# rm(df_beta, map)
# gc()


#---- Read correlation matrix previously computed ----
# tmp <- paste(sub(".rds", "", map_ld_rds), pheno, sep="_")
# tmp <- basename(tmp)
# tmp <- file.path(tmpdir, tmp)
tmp <- tempfile(tmpdir = tmpdir)
corr <- NULL
for (mychr in 1:22){
  ix <- df_beta_good[chr == mychr, `_NUM_ID_`]
  
  # retrieve indexes relative to the correlation matrix positions
  ix2 <- match(ix, which(mapLDref$chr == mychr))
  
  # Get correlation file
  fi <- grepl(glue("chr{mychr}.rds"), corfiles)
  ff <- corfiles[fi]
  corr_chr <- readRDS(ff)[ix2, ix2]
  
  if (is.null(corr)) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# Remove the following line for local execution
# tmpdir <- "tmp-data"
# 
# The following block has been substituted by the script
# "bin/create_cor.R
# >>>>>> BEGIN BLOCK -------------------------------------
# ldrefdir <- dirname(ldfiles[1])
# tmp <- tempfile(tmpdir = tmpdir)
# t1 <- Sys.time()
# for (chr in 1:22){
#   tc1 <- Sys.time()
#   cat(chr, ".. ", sep = "")
#   
#   ## indices in 'df_beta'
#   ind.chr <- which(df_beta_good$chr == chr)
#   ## indices in 'map_ldref'
#   cat("Ind.chr2...\n")
#   ind.chr2 <- df_beta_good$`_NUM_ID_`[ind.chr]
#   ## indices in 'corr_chr'
#   cat("Ind_chr3...\n")
#   ind.chr3 <- match(ind.chr2, which(map$chr == chr))
#   
#   #---- Read correlation matrix ----
#   idf <- match(glue("{ldrefdir}/ld_chr{chr}.rds"), ldfiles)
#   if (!is.na(idf)){
#     cat("read matrix...\n")
#     corr_chr <- readRDS(ldfiles[[idf]])[ind.chr3, ind.chr3]
#     
#     cat("Compute ld by snp..\n")
#     tmpld <- Matrix::colSums(corr_chr^2)
#     df_beta_good[ind.chr[!is.na(ind.chr3)], ld:=tmpld]
#     
#     if (chr == 1) {
#       corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
#     } else {
#       corr$add_columns(corr_chr, nrow(corr))
#     }
#   } else {
#     stop(glue("No such file or directory: .... ld_ref/ld_chr{chr}.rds"))
#   }
#   tc2 <- Sys.time()
#   cat(glue("Chromosome {chr} computed in: {format(tc2-tc1)}"), "\n")
# }
# t2 <- Sys.time()
# cat(glue("Compute correlation in {format(t2-t1)}...\n"))
# <<<<<< END BLOCK -------------------------------------

#---- Export df_beta_good for future computation ----
saveRDS(df_beta_good, file=df_beta_good_rds)

#---- Heritability estimation of LD score regression ----
# to be used as a starting value in LDpred2-auto
t1 <- Sys.time()
cat("Running heritability estimation...")
# (ldsc <- with(df_beta_good, snp_ldsc(ld, ld_size = nrow(map),
#                                 chi2 = (beta / beta_se)^2,
#                                 sample_size = n_eff,
#                                 ncores = NCORES, blocks = NULL)))
(ldsc <- snp_ldsc(df_beta_good$ld, nrow(df_beta_good), 
                  chi2 = (df_beta_good$beta / df_beta_good$beta_se)**2,
                  ncores = NCORES, 
                  sample_size = mean(df_beta_good$n_eff, na.rm=TRUE))
)

t2 <- Sys.time()
# This estimation is to high and gives instability of the models, so we 
# try with a different (manual) value
h2_est <- ldsc[["h2"]]
cat(glue("run in {format(t2-t1)}"), "\n")
cat(glue("Estimated heritability: {h2_est}"), "\n")
cat(glue("Saving {heritability_file}..."), "\n")
saveRDS(ldsc, file=heritability_file)


# LDpred2-auto
# This model always returns NAs...
# So let's try with grid-model
if (file.exists(model_file)){
  multi_auto <- readRDS(model_file)
} else {
  t1 <- Sys.time()
  cat("Running LDpred-auto estimation...\n")
  multi_auto <- snp_ldpred2_auto(corr, df_beta_good, h2_init = h2_est,
                                 vec_p_init = seq_log(1e-4, 0.9, length.out = 10),
                                 allow_jump_sign = FALSE, shrink_corr = 0.95,
                                 ncores = NCORES) # 5 min
  t2 <- Sys.time()
  cat(glue("LDpred model run in {format(t2-t1)}"),"\n")
  cat(glue("Saving {model_file}..."), "\n")
  saveRDS(multi_auto, file=model_file)
}

range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

if (length(keep) == 0){
  keep <- 1:length(multi_auto)
}
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

# Save beta for the LDPred auto model
saveRDS(beta_auto, file=beta_file)

#---- LDpred2 - grid ----
cat("Running ldpred2-grid model\n")
# (h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4))
(h2_seq <- round(h2_est * c(0.01, 0.3, 0.7, 1, 1.4, 2), 4))
h2_seq <- h2_seq[h2_seq > 0]
(p_seq <- signif(seq_log(1e-5, 1, length.out = 10), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))
saveRDS(params, file=params_grid_file)

beta_grid <- snp_ldpred2_grid(corr, df_beta_good, params, ncores = NCORES)#, ind.corr = df_beta_good$`_NUM_ID_`)
saveRDS(beta_grid, file=beta_grid_file)

# TODO: move the prediction into another script
pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta_good[["_NUM_ID_.ss"]])
saveRDS(pred_grid, file=pred_grid_file)

# Quit
quit(save="no")

#---- Lassosum2 computation ----
# This has been moved into another script!
cat("Computing beta with lassosum2...")
t1 <- Sys.time()
beta_lassosum2 <- snp_lassosum2(corr, df_beta_good, ncores = NCORES)
t2 <- Sys.time()
cat(glue("in {format(t2 - t1)}"), "\n Saving...\n")
saveRDS(beta_lassosum2, file=beta_file)

#---- Lassosum2 parameters ----
(params2 <- attr(beta_lassosum2, "grid_param"))
saveRDS(params2, file=params_grid_file)



cat("Computing best beta...\n")
# Filter for best chains and average remaining ones
# -> the effects sizes of your polygenic score
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- (range > (0.95 * quantile(range, 0.95))))
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
cat(glue("Saving {beta_file}...\n"))
saveRDS(beta_auto, file=beta_file)