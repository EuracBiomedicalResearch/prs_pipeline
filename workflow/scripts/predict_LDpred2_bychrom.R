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

#---- Input/Output ----
# TODO: add data.frame of input variants used for the computation
beta_auto_file <- snakemake@input[["pred_auto"]]

pred_file <- snakemake@output[["pred_rds"]]
pred_csv <- snakemake@output[["pred_csv"]]
map_out <- snakemake@output[["map_file"]]

#---- Resources ----
NCORES <- snakemake@threads

#---- Convert to numeric ----
# This function converts to numeric and set NA to 0
to_double <- function(x){
    y <- as.numeric(x)
    y[is.na(y)] <- 0
    return(y)
}

#---- Cycle over all files and sum prs up ----
i <- 1
for (rr in beta_auto_file){
    df <- readRDS(rr)
    if (i == 1){
        prs <- to_double(df$PRS)
    } else {
      prs <- prs + to_double(df$PRS)
    }
    i <- i + 1
}
pred_auto <- df[, .SD, .SDcols =!c("PRS")]
pred_auto$PRS <- prs

#---- Save predictions into a rdsfile ----
cat("Saving prediction to file rds...\n")
saveRDS(pred_auto, file = pred_file)

#---- Save predictions into a csv file ----
cat("Saving prediction to file csv...\n")
write.table(pred_auto, file = pred_csv, sep = "\t",
    row.names = FALSE, col.names = TRUE)
