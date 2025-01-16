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

# Output
pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
params_csv <- snakemake@output[["params_csv"]]

for (i in 1:length(predfiles)){
    prs <- readRDS(predfiles[i])
    if (i == 1){
        tmpprs <- prs[, 3:ncol(prs)]
    } else {
        tmpprs <- tmpprs + prs[, 3:ncol(prs)]
    }
}

# Add family ID and sample ID
prs_all_samps <- cbind(prs[, 1:2], tmpprs)

# Save PRS RDS
saveRDS(prs_all_samps, file=pred_file)

# Save PRS csv
write.table(prs_all_samps, file=pred_csv, col.names = TRUE, row.names = FALSE, 
sep = "\t")

# Read grid parameters
chr <- sub(".+chr(\\d+).+", "\\1", basename(predfiles[i]))
clump_res_file <- file.path(dirname(predfiles[i]), 
    glue("multi_prs_ct_chr{chr}.rds"))
clump_res <- readRDS(clump_res_file)

grid_params <- attr(attr(clump_res, "all_keep"), "grid") %>%
  mutate(thr.lp = list(attr(clump_res, "grid.lpS.thr"))) %>%
  unnest(cols = "thr.lp") %>% 
  mutate(col_id = 1:n())

# Save parameters to a csv file
write.table(grid_params, file=params_csv, col.names = TRUE, row.names = FALSE, 
sep = "\t")