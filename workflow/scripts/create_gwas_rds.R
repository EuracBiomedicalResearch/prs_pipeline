#!/shared/bioinf/R/bin/Rscript-4.3-BioC3.17
# Create RDS for Snakemake pipeline
library(data.table)
library(jsonlite)

#---- Get parameter from SNAKEMAKE ----
infile <- snakemake@input[["gwas_file"]]
outfile <- snakemake@output[["gwas_rds"]]
col_dict <- snakemake@params[["col_dict"]]
print(col_dict)
# colnam <- jsonlite::parse_json(col_dict)


#---- Read files ----
gwas_data <- fread(infile)
setnames(gwas_data, names(col_dict), as.character(col_dict))

# gwas_data <- setNames(gwas_data, c("chr", "pos", "rsid", "a0", "a1", "freq", "beta", "beta_se", "p", "n_eff"))

#---- Save RDS ----
saveRDS(gwas_data, file=outfile)
