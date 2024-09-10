#!/shared/bioinf/R/bin/Rscript-4.3-BioC3.17
# Create RDS for Snakemake pipeline
library(data.table)
library(jsonlite)
library(glue)

#---- Get parameter from SNAKEMAKE ----
infile <- snakemake@input[["gwas_file"]]
outfile <- snakemake@output[["gwas_rds"]]
gwas_conf <- snakemake@params[["gwas_conf"]]
res_dir <- snakemake@params[["res_dir"]]
col_dict <- gwas_conf[["format"]]

# Read name translation from `GWASlab/formatbooks` to `bigsnpr`
fmt_2_R <- jsonlite::read_json(file.path(res_dir, "formatbook_to_R.json"))$format

#---- Check if names of the columns are available ----
if (gwas_conf[["header"]] == TRUE){
  if (is.null(col_dict)){
    # Throw error because no column names/format is provided
    stop(glue("Please provide a column schema for GWAS: {infile}"))
  } else {
      if (is.list(col_dict)){
        # Column definition provided by the user
        format_names <- col_dict[["format_dict"]]
    } else if (is.character(col_dict) & col_dict != ""){
        # Column definition provided by formatbook through predefined formats
        format_names <- jsonlite::read_json(file.path(res_dir, "formatbook", "formatbook.json"))
        format_names <- format_names[[col_dict]][["format_dict"]]
    } else {
      stop(glue("Please provide a column schema for GWAS: {infile}"))
    }
  }
}

#---- Read GWAS sumstat file ----
gwas_data <- fread(infile)

#---- Handle cases when passing format for column names ----
if (is.character(col_dict)){ # & col_dict == "auto"){
  ix <- match(names(gwas_data), names(format_names))
  ix <- ix[!is.na(ix)]
  format_names <- as.list(as.data.frame(format_names)[,ix])
}

#---- Set names ----
# 1. Names to GWAS lab standard
if (!is.null(format_names)){
  setnames(gwas_data, names(format_names), as.character(format_names), skip_absent = TRUE)
}
# 2. From GWASlab to bigsnpr
setnames(gwas_data, names(fmt_2_R), as.character(fmt_2_R), skip_absent = TRUE)

#---- Change chromosome name if genome_build is hg38 ----
# Column: chr
#         "chr1" -> "1"
# Plink do not use chr in front of the chromosome code
if (gwas_conf$genome_build == "hg38"){
  gwas_data[, chr:=as.character(chr)]
  ck_chrom <- grepl("chr\\d{1,2}", gwas_data$chr, perl=TRUE)
  gwas_data[ck_chrom, chr:=sub("chr", "", chr)]
  gwas_data[, chr:=as.integer(chr)]
}

# TODO: Add filter for Minor allele frequency

# TODO: Compute column P of p-values when only log10p is provided and the opposite

# TODO: Add n_eff when binary phenotype is available and N_Cases and N_Controls are available
# dirty solution provided here based on https://privefl.github.io/bigsnpr/articles/LDpred2.html,
# something better should be provided...
if (!any(names(gwas_data) == "n_eff")){
  n <- quantile(8 / gwas_data$beta_se^2, 0.999)
  gwas_data[, n_eff:=n]
}

#---- Save RDS ----
saveRDS(gwas_data, file=outfile)
