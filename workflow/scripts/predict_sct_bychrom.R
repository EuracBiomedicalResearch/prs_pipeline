#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr)
library(rmio)
library(tidyr)

#---- Setup ----
# Input
predfiles <- snakemake@input[["pred_chrom"]]
paramsfiles <- snakemake@input[["params_chrom"]]
info_snp_files <- snakemake@input[["info_snp"]]

# Output
pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
params_csv <- snakemake@output[["params_csv"]]
info_snp_csv <- snakemake@output[["info_snp_csv"]]

#---- Create parameter data.frame ----
dfannot <- c()
j <- 0
for (i in 1:length(paramsfiles)){
  f <- paramsfiles[i]
  if (file.size(f) > 0){
    j <- j + 1
    df <- read.csv(f, header = TRUE, sep = "\t")
  }
  if (j == 1) {
    dfannot <- df
  }
}
nparams <- nrow(dfannot)

#---- Create data.frame with info on variants used ----
info_snp_list <- list()
for (i in 1:length(info_snp_files)){
  f <- info_snp_files[i]
  if (file.size(f) > 0){
    info_snp_list[[i]] <- fread(f, header=TRUE)
  }
}
info_snp <- data.table::rbindlist(info_snp_list)

#---- Collect prediuctions ----
for (i in 1:length(predfiles)){
  prs <- readRDS(predfiles[i])
  prsvals <- prs[, 3:ncol(prs)]

  # Handle case when PRS is null
  # Should do this way in case the PRS is null for the first chromosome....
  if (ncol(prsvals) < nparams){
    nn <- paste("p", 1:nparams)
    prsvals <- data.frame(sapply(nn, function(x){rep(0, nrow(prs))}))
  }
  if (i == 1){
    tmpprs <- prsvals
  } else {
    tmpprs <- tmpprs + prsvals
  }
}

# Add family ID and sample ID
prs_all_samps <- cbind(prs[, 1:2], tmpprs)
names(prs_all_samps)[1:2] <- c("family.ID", "sample.ID")

#---- Saving outputs ----
# Save PRS RDS
saveRDS(prs_all_samps, file=pred_file)

# Save PRS csv
write.table(prs_all_samps, file=pred_csv, col.names = TRUE, row.names = FALSE, 
sep = "\t")

# Save info_snp
write.table(info_snp, file=info_snp_csv, row.names = FALSE, col.names = TRUE, 
  sep = "\t", quote=FALSE)

# Save parameters
write.table(dfannot, file=params_csv, row.names = FALSE, col.names = FALSE,
  sep="\t")