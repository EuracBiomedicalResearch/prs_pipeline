#!/shared/bioinf/R/bin/Rscript-4.3-BioC3.17
#
#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)#, lib.loc = "/home/mfilosi/R/rocker-rstudio/4.2")

Sys.setenv(OPENBLAS_NUM_THREADS=1)

#---- Setup ----
bed <- snakemake@input[["bed"]]
bim <- snakemake@input[["bim"]]
fam <- snakemake@input[["fam"]]

ofile <- snakemake@output[[1]]
ofilebk <- sub(".rds$", "", ofile)
mapfile <- snakemake@output[[3]]

# Resources
ncores <- snakemake@threads
tmpdir <- snakemake@resources[["tmpdir"]]

#---- Read Genotype using bigsnpr ----
cat(glue("Read bed file {bed} using backingfile {ofilebk}\n"))
t1 <- Sys.time()
geno <- snp_readBed(bed, backingfile = ofilebk)
t2 <- Sys.time()
cat(glue("Read in {format(t2-t1)}\n"), "\n")

#---- Attach SNP info ----
cat(glue("Attach genotype from file {ofile}\n"), "\n")
geno <- snp_attach(ofile)
t1 <- Sys.time()
cat(glue("Attached in {format(t1 - t2)}\n"), "\n")
map <- geno$map
map <- setNames(map[-3], c("chr", "rsid", "pos", "a1", "a0"))
CHR <- as.integer(map$chr)
POS <- map$pos

#---- Genetic Position ----
# options(bigstatsr.check.parallel.blas = TRUE)
# bigparallelr::set_blas_ncores(1)
# 
# cat("Compute genetic positions...\n")
# t1 <- Sys.time()
# POS2 <- snp_asGeneticPos(CHR, POS, dir = tmpdir, ncores = ncores)
# map$genetic_pos <- POS2
# t2 <- Sys.time()
# cat(glue("Genetic distance computed in {format(t2-t1)}\n"))

#---- Compute Allele frequency ----
cat("Compute variant frequency...\n")
freq <- big_colstats(geno$genotypes, ncores = ncores)$sum / (2 * nrow(geno$genotypes))
map$freq <- freq
t1 <- Sys.time()
cat(glue("Variant frequency computed in {format(t1 -t2)}"), "\n")

#---- Re-read the plink files with selected genotypes ----
if (file.exists(glue("{ofilebk}.bk"))){
  cat("Removing file bk..\n`")
  file.remove(glue("{ofilebk}.bk"))
}
ix <- which(freq > 0.01 & freq < 0.99)
geno <- snp_readBed2(bed, backingfile = ofilebk, ind.col = ix, ncores=ncores)
cat(glue("Read from the plink files second time.."), "\n")

#---- Saving map files after filtering ----
cat("Saving map file after filtering...\n")
saveRDS(map[ix, ], file=mapfile)

#---- Update Genotypes ----
# cat("Select only variants with freq higher than 0.01\n")
# ix <- which(freq > 0.01 & freq < 0.99)
# cat("copying the matrix...")
# t1 <- Sys.time()
# G2 <- big_copy(geno$genotypes, ind.col=ix)
# t2 <- Sys.time()
# cat(glue(" in {format(t2 - t1)}"), "\n")
# geno$genotypes <- G2
# geno$map <- map[ix, ]
# 
# map <- map[ix, ]
# 
# cat("Saving genotype file after filtering")
# snp_save(geno)

