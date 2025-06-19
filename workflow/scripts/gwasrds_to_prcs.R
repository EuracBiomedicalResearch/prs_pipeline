#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(data.table)

#---- Input/Output ----
gwas_rds <- snakemake@input[["gwas_rds"]]
bim_file <- snakemake@input[["bim_file"]]
chrom <- snakemake@wildcards[["chrom"]]


#---- Read in the GWAS ----
cat("Reading file...\n")
gwas <- readRDS(gwas_rds)
gwas <- gwas[chr == chrom, ]


if (sum(names(gwas) == "RSID") == 0){
    #---- Read bim ----
    bim <- data.table::fread(bim_file)
    names(bim) <- c("chr", "BIMID", "cM", "pos", "a0", "a1")

    # Match snps
    info_snp <- snp_match(gwas, bim, match.min.prop=0)  
    gwas <- info_snp[, c("BIMID", "a0", "a1", "beta", "beta_se", "freq")]
    names(gwas)[1] <- "snpid"
} else {
    gwas[, snpid:=RSID]
}

#---- Extract only needed columns ----
gwas <- gwas %>% 
    dplyr::filter(freq > 0.01, freq < 0.99) %>% 
    dplyr::select(any_of(c("snpid", "a0", "a1", "beta", "beta_se")))
names(gwas) <- c("SNP", "A1", "A2", "BETA", "SE")


#---- Write gwas output ----
cat("Writing file...\n")
data.table::fwrite(gwas, file=snakemake@output[["gwas_prscs"]], sep="\t")
