
rule prs_CT:
  message:
    "Running SCT PRS"
  input:
    geno_data = os.path.join(geno_dir, "qc_geno_all.rds"),
    map_data = os.path.join(geno_dir, "qc_geno_all_map.rds"),
    gwas_data= os.path.join(odir, "gwas.rds"),
    # train_data = "data/train_pheno_file.rds"
  threads: 16
  resources:
    mem_mb = 64000
  output:
    clump_opt = os.path.join(odir, "sct/clump_res_ct.rds"),
    multi_prs = os.path.join(odir, "sct/multi_prs_ct.rds"),
    multi_prs_bk = os.path.join(odir, "sct/multi_prs_ct.bk"),
    pred_rds = os.path.join(odir, "sct/prs.rds"),
    pred_csv = os.path.join(odir, "sct/prs.csv")
  params:
      force = 'FALSE'
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_PRS_CT_v2.R"
