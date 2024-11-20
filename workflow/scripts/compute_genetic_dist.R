#!/shared/bioinf/R/bin/Rscript-4.3-BioC3.17
#
#---- Libraries ----
library(dplyr)
library(glue)
library(bigsnpr)

Sys.setenv(OPENBLAS_NUM_THREADS=1)

#---- Setup ----
mapfile <- snakemake@input[["mapfile"]]

#---- Get parameters ----
tmpdir <- snakemake@resources[["tmpdir"]]
ncores <- snakemake@threads

#---- Set BLAS parallel env to 1 ----
options(bigstatsr.check.parallel.blas = TRUE)
bigparallelr::set_blas_ncores(1)

#---- Read Genotype using bigsnpr ----
cat("Reading map file...\n")
map <- readRDS(mapfile)
CHR <- as.integer(map$chr)
POS <- map$pos

#---- Compute genetic position ----
cat("Compute genetic position...")
t1 <- Sys.time()
POS2 <- snp_asGeneticPos(CHR, POS, dir = tmpdir, ncores=ncores)
t2 <- Sys.time()
map$genetic_pos <- POS2
cat(glue(" in {format(t2 - t1)}"), "\n")

#---- Compute Allele frequency ----
# freq <- big_colstats(geno$genotypes, ncores = ncores)$sum / (2 * nrow(geno$genotypes))
# map$freq <- freq

#---- Save output ----
cat("Saving output...")
saveRDS(map, file=snakemake@output[[1]])
