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

# Output
pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
params_csv <- snakemake@output[["params_csv"]]


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

for (i in 1:length(predfiles)){
  prs <- readRDS(predfiles[i])
  prsvals <- prs[, 3:ncol(prs)]

  # Handle case when PRS is null
  # Should do this way in case the PRS is null for the first chromosome....
  if (ncol(prsvals) < nparams){
    nn <- paste("V", 1:nparams)
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

# Save PRS RDS
saveRDS(prs_all_samps, file=pred_file)

# Save PRS csv
write.table(prs_all_samps, file=pred_csv, col.names = TRUE, row.names = FALSE, 
sep = "\t")


# This is a temporary line just to let the work finish correctly
# before solving the bug describe in the TODO
file.create(params_csv)

# TODO: improve this part of code giving the params grid as output
# # Read grid parameters
# chr <- sub(".+chr(\\d+).+", "\\1", basename(predfiles[i]))
# clump_res_file <- file.path(dirname(predfiles[i]), 
#     glue("multi_prs_ct_chr{chr}.rds"))
# clump_res <- readRDS(clump_res_file)

# grid_params <- attr(attr(clump_res, "all_keep"), "grid") %>%
#   mutate(thr.lp = list(attr(clump_res, "grid.lpS.thr"))) %>%
#   unnest(cols = "thr.lp") %>% 
#   mutate(col_id = 1:n())

# # Save parameters to a csv file
# write.table(grid_params, file=params_csv, col.names = TRUE, row.names = FALSE, 
# sep = "\t")