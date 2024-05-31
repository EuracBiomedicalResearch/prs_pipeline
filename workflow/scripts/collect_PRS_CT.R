#---- Libraries ----
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(glue)
library(bigsnpr, lib.loc = "/home/mfilosi/R/rocker-rstudio/4.2")
library(rmio)


multiprs_files <- snakemake@input[["multiprs_file"]]

for (f in multiprs_files){
  df <- readRDS(f)
}