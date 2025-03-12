#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(data.table)
library(tidyr)

#---- Input/output ----
prs_pred_file <- snakemake@input[["pred_rds"]]

prs_dist_file <- snakemake@output[["prs_dist"]]

#---- Read RDS ----
prs_df <- readRDS(prs_pred_file)

#---- Create plot ----
dist_plot <- prs_df %>% pivot_longer(c(-family.ID, -sample.ID)) %>% 
    ggplot(aes(x=value, group=name, color=name)) + 
    geom_density(show.legend = FALSE) + theme_minimal() + 
    xlab("PRS distribution")

#---- Saving plot ----
ggsave(dist_plot, file=prs_dist_file)
