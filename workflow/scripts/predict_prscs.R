#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)#, lib.loc = "/home/mfilosi/R/rocker-rstudio/4.2")
library(data.table)

# Avoid nested parallel computation...
Sys.setenv(OPENBLAS_NUM_THREADS=1)

#---- Input/Output ----
# Input
bim_files <- snakemake@input[["bim"]]
beta_file <- snakemake@input[["beta_shrinked"]]
geno_file <- snakemake@input[["genotype_rds"]]
# Output
map_file <- snakemake@output[["map_file"]]
pred_file <- snakemake@output[["pred_file"]]
pred_csv <- snakemake@output[["pred_csv"]]
# Params
ld_dir <- snakemake@params[["ld_dir"]]
lddata <- snakemake@params[["lddata"]]
gwas_conf <- snakemake@params[["gwas_conf"]]

#---- Read genome build for the GWAS ----
gwas_build = gwas_conf[["genome_build"]]

# -----------------------------------------------------
# Uncomment the following lines if 
# you want to run in local version not with snakemake
# -----------------------------------------------------
# bim_files <- sapply(1:22, function(x){glue("data/geno/chris/qc_geno_chr{x}_rsid.bim")})
# pred_file <- "results/prscs/pred_prscs.rds"
# map_file <- "results/prscs/map_prscs.rds"
# beta_file <- "results/prscs/egfr_pst_all.txt"

#---- Read in bim files ----
bimdflist <- list()
for (f in bim_files){
  bimdflist[[f]] <- data.table::fread(f, col.names = c("chr", "id", "cm", "pos", "a0", "a1"))
  # bimdf <- rbind(bimdf, tmpdf)
}
bimdf <- data.table::rbindlist(bimdflist)

# Order by chromosome and position
setorder(bimdf, chr, pos)

# Create original SNPID 
# bimdf[, ID:=paste(CHROM, BP, A1, A0, sep=":")]

#---- Read in beta file ----
betas <- data.table::fread(beta_file, header=FALSE, 
                           col.names=c("chr", "id", "pos", "a0", "a1", "beta"))

# If build is 38 need to adjust position according to the SNPinfo file
if (gwas_build == "hg38"){
  snpinfo <- data.table::fread(file.path(ld_dir, paste0("snpinfo_", lddata, "_hm3_hg38")))
  snpinfo <- snpinfo[, -c("CHR", "A1", "A2", "MAF")]
  betas <- merge(betas, snpinfo, by.x="id", by.y="SNP", sort=FALSE)
  betas <- betas[, -"pos"]
  setnames(betas, "BP", "pos")
}

# bimdf[, betas:=NA]
# bimdf_betas <- merge(bimdf, betas, by=c("CHROM", "RSID", "BP", "A0", "A1"))
# setorder(bimdf, chr, pos)

#---- Attach genotype file ----
# geno_file <- "data/geno/chris/qc_geno_all.rds"
geno <- snp_attach(geno_file)

#---- Mapping geno file with new annotated bim ----
betas_match <- snp_match(sumstats = betas, info_snp = bimdf)
# idx <- match(geno$map$marker.ID, bimdf_betas$ID)
# genocol <- which(!is.na(idx))

prscs_pred <- big_prodVec(geno$genotypes,
                          betas_match[["beta"]], 
                          ind.col = betas_match[["_NUM_ID_"]])

#---- Extract map beta after shrinkage and correction ----
# df_beta <- geno$map
# df_beta$numid <- 1:nrow(df_beta)
# df_beta <- df_beta[genocol,]
# df_beta$beta <- bimdf_betas[["BETA"]][idx[!is.na(idx)]]
# df_beta$rsid2 <-  bimdf_betas[["RSID"]][idx[!is.na(idx)]]

saveRDS(betas_match, file=map_file)

#---- Save predictions ----
# exporttdf <- data.table(AID=geno$fam$sample.ID, prs_prscs=prscs_pred)
prs <- data.table(family.ID=geno$fam$family.ID, sampleID=geno$fam$sample.ID, prs_prscs=prscs_pred)
# save RDS
saveRDS(prs, file=pred_file)

# Save CSV
write.table(prs, file=pred_csv, col.names=TRUE, row.names=FALSE, sep="\t")