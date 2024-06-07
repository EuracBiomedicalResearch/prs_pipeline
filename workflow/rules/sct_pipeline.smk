
rule prs_CT:
  message:
    "Running SCT PRS"
  input:
    geno_data = "data/geno/qc_geno_all.rds",
    map_data = "data/geno/qc_geno_all_map.rds",
    gwas_data= "results/{pheno}/gwas.rds",
    # train_data = "data/train_pheno_file.rds"
  threads: 16
  resources:
    mem_mb = 64000
  output:
    clump_opt = "results/{pheno}/sct/clump_res_ct.rds",
    multi_prs = "results/{pheno}/sct/multi_prs_ct.rds",
    multi_prs_bk = "results/{pheno}/sct/multi_prs_ct.bk",
    pred_rds = "results/{pheno}/sct/prs.rds",
    pred_csv = "results/{pheno}/sct/prs.csv"
    # final_mod = "results/{pheno}/final_mod_CT.rds",
    # cov_pred = 'results/covariate_predictions.rds',
    # gen_pred = 'results/genetic_predictions.rds',
    # pred = 'results/predictions.rds'
  params:
      force = 'FALSE'
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_PRS_CT_v2.R"
