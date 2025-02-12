rule estimate_heritability:
  message:
    "Estimate heritability genome wide"
  input:
    map_ld_rds = ancient(hm3map),
    gwas_rds = os.path.join(odir, "gwas.rds")
  output:
    heritability_file = os.path.join(odir, "ldpred2/heritability_ldpred2.rds")
  params:
    genotype_conf = genotype_conf
  resources:
    mem_mb = get_mem_mb
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/estimate_heritability.R"

rule prs_LDpred2:
  input:
    # train_data = 'data/train_pheno_file.rds',
    map_rds = os.path.join(geno_dir, "qc_geno_chr{chrom}_map.rds"),
    gwas_rds = os.path.join(odir, "gwas.rds"),
    genotype_rds = os.path.join(geno_dir, "qc_geno_chr{chrom}.rds"),
    heritability = os.path.join(odir, "ldpred2/heritability_ldpred2.rds"),
    map_ld_rds = ancient(hm3map),
    corfiles = ancient(os.path.join(hm3path, "ldref_hm3_plus", "LD_with_blocks_chr{chrom}.rds"))
    # corfiles = ancient(hm3corr)
    # map_ld_rds = "resources/ld_ref/map_hm3_plus.rds",
    # corfiles = expand("resources/ld_ref/ldref_hm3_plus/LD_with_blocks_chr{chrom}.rds", 
    #                  chrom=range(1,23))
  output:
    # final_mod = os.path.join(odir, "ldpred2/final_mod_chr{chrom}_ldpred2.rds"),
    # beta_auto = os.path.join(odir, "ldpred2/beta_auto_chr{chrom}_ldpred2.rds"),
    df_good = os.path.join(odir, "ldpred2/df_beta_good_with_ld_chr{chrom}.rds"),
    # beta_grid = os.path.join(odir, "ldpred2/beta_grid_chr{chrom}_ldpred2.rds"),
    pred_grid = os.path.join(odir, "ldpred2/pred_grid_chr{chrom}_ldpred2.rds"),
    params_grid = os.path.join(odir, "ldpred2/params_grid_chr{chrom}_ldpred2.rds"),
    pred_auto = os.path.join(odir, "ldpred2/pred_auto_chr{chrom}_ldpred2.rds")
  resources:
    mem_mb = get_mem_mb,
    tmpdir = "tmp-data"
  params:
    genotype_conf = genotype_conf
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_PRS_ldpred2_bychrom.R"


rule predict_LDpred2:
  input:
    # train_data = 'data/train_pheno_file.rds',
    # map_rds = 'data/geno/rds/qc_geno_all_map.rds',
    # map_rds = os.path.join(geno_dir, "qc_geno_all_map.rds"),
    # gwas_rds = os.path.join(odir, "gwas.rds"),
    # map_good = os.path.join(odir, "ldpred2/df_beta_good_with_ld.rds"),
    # genotype_rds = os.path.join(geno_dir, "qc_geno_all.rds"),
    # beta_auto = os.path.join(odir, "ldpred2/beta_auto_ldpred2.rds")
    pred_auto = expand_chrom(os.path.join(odir, "ldpred2/pred_auto_chr{chrom}_ldpred2.rds"))
  output:
    pred_rds = os.path.join(odir, "ldpred2/prs.rds"),
    pred_csv = os.path.join(odir, "ldpred2/prs.csv"),
    # map_file = os.path.join(odir, "ldpred2/map_prs.rds")
  resources:
    mem_mb = get_mem_mb
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/predict_LDpred2_bychrom.R"

rule prs_lassosum2:
  message:
    "Run PRS estimation with lassosum2"
  input:
    map_rds = os.path.join(geno_dir, "qc_geno_chr{chrom}_map.rds"),
    gwas_rds = os.path.join(odir, "gwas.rds"),
    genotype_rds = os.path.join(geno_dir, "qc_geno_chr{chrom}.rds"),
    map_ld_rds = ancient(hm3map),
    corfiles = ancient(os.path.join(hm3path, "ldref_hm3_plus", "LD_with_blocks_chr{chrom}.rds"))
    # map_ld_rds = "resources/ld_ref/map_hm3_plus.rds", 
    # corfiles = expand("resources/ld_ref/ldref_hm3_plus/LD_with_blocks_chr{chrom}.rds", 
                     # chrom=range(1,23))
  output:
    # beta_file = os.path.join(odir, "lassosum2/beta_lassosum2_chr{chrom}.rds"),
    df_beta = os.path.join(odir, "lassosum2/df_beta_good_lassosum2_chr{chrom}.rds"),
    params_file = os.path.join(odir, "lassosum2/params_grid_lassosum2._chr{chrom}rds"),
    pred_file = os.path.join(odir, "lassosum2/pred_lassosum2_chr{chrom}.rds")
  params:
    genotype_conf = genotype_conf
  resources:
    mem_mb = get_mem_mb,
    tmpdir = "tmp-data"
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_lassosum2_bychrom.R"


rule predict_lassosum2:
  message:
    "Predict PRS with lassosum2"
  input:
    # genotype_rds = os.path.join(geno_dir, "qc_geno_all.rds"),
    # df_beta = os.path.join(odir, "lassosum2/df_beta_good_lassosum.rds"),
    # beta_lassosum = os.path.join(odir, "lassosum2/beta_lassosum2.rds"),
    pred_lassosum = expand_chrom(os.path.join(odir, "lassosum2/pred_lassosum2_chr{chrom}.rds"))
  output:
    # pred = "results/{pheno}/lassosum2/pred_lassosum2.rds"
    pred_rds = os.path.join(odir, "lassosum2/prs.rds"),
    pred_csv = os.path.join(odir, "lassosum2/prs.csv"),
    # map_file = os.path.join(odir, "lassosum2/map_prs.rds"),
  resources:
    mem_mb=get_mem_mb
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/predict_lassosum2_bychrom.R"