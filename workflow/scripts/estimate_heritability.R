#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(data.table)

#---- Input/Ouput ----
# Input
gwas_rds <- snakemake@input[["gwas_rds"]]
map_ld_rds <- snakemake@input[["map_ld_rds"]]
genotype_conf <- snakemake@params[["genotype_conf"]]
# Output
heritability_file <- snakemake@output[["heritability_file"]]
genotype_conf <- snakemake@params[["genotype_conf"]]

#---- Resources ----
# Avoid nested parallel computation...
Sys.setenv(OPENBLAS_NUM_THREADS=1)
NCORES <- snakemake@threads

#---- Load gwas ----
gwas <- readRDS(gwas_rds)

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
df_beta_good <- snp_match(gwas, mapLDref, match.min.prop = 0)
df_beta_good <- as.data.table(df_beta_good)

#---- Heritability estimation of LD score regression ----
# to be used as a starting value in LDpred2-auto
t1 <- Sys.time()
cat("Running heritability estimation...\n")
(ldsc <- snp_ldsc(df_beta_good$ld, nrow(df_beta_good), 
                  chi2 = (df_beta_good$beta / df_beta_good$beta_se)**2,
                  ncores = NCORES, 
                  sample_size = mean(df_beta_good$n_eff, na.rm=TRUE))
)
t2 <- Sys.time()
cat(glue("run in {format(t2-t1)}"), "\n")
cat(glue("Estimated heritability: {ldsc[['h2']]}"), "\n")
cat(glue("Saving {heritability_file}..."), "\n")
saveRDS(ldsc, file=heritability_file)