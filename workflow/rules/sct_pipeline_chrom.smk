rule prs_CT:
  message:
    "Running SCT PRS"
  input:
    geno_data = os.path.join(geno_dir, "qc_geno_chr{chrom}.rds"),
    map_data = os.path.join(geno_dir, "qc_geno_chr{chrom}_map.rds"),
    gwas_data = os.path.join(odir, "gwas.rds")
  resources:
    mem_mb = get_mem_mb
  output:
    # clump_opt = os.path.join(odir, "sct/clump_res_ct_chr{chrom}.rds"),
    # multi_prs = os.path.join(odir, "sct/multi_prs_ct_chr{chrom}.rds"),
    # multi_prs_bk = os.path.join(odir, "sct/multi_prs_ct_chr{chrom}.bk"),
    pred_rds = os.path.join(odir, "sct/prs_chr{chrom}.rds"),
    pred_csv = os.path.join(odir, "sct/prs_chr{chrom}.csv"),
    params_csv = os.path.join(odir, "sct/params_chr{chrom}.csv")
  params:
      force = "FALSE"
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_PRS_CT_v2.R"

rule predict_CT:
  message:
    "Running prediction using CT algorithm"
  input:
    pred_chrom = expand_chrom(os.path.join(odir, "sct/prs_chr{chrom}.rds")),
    params_chrom = expand_chrom(os.path.join(odir, "sct/params_chr{chrom}.csv"))
  output:
    pred_rds = os.path.join(odir, "sct/prs.rds"),
    pred_csv = os.path.join(odir, "sct/prs.csv"),
    params_csv = os.path.join(odir, "sct/prs_params.csv")
  resources:
    mem_mb = get_mem_mb
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/predict_sct_bychrom.R"