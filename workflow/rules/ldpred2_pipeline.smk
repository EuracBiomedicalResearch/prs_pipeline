# include: "preprocess.smkg"

# rule allldref:
#   input:
#     expand('results/{dataset}/heritability_ldpred2.rds', dataset=['hm3', 'chris']),
#     expand('results/{dataset}/final_mod_ldpred2.rds', dataset=['hm3', 'chris']),
#     expand('results/{dataset}/beta_auto_ldpred2.rds', dataset=['hm3', 'chris']),
#     expand('data/{dataset}/df_beta_good_with_ld.rds', dataset=['hm3', 'chris']),
#     expand('results/{dataset}/beta_grid_ldpred2.rds', dataset=['hm3', 'chris'])
#     # expand('results/{dataset}/pred_grid_ldpred2.rds', dataset=['chris', 'hm3'])

# rule all_lasso:
#   input:
#     expand('results/{dataset}/beta_lassosum2.rds', dataset=[ 'chris', 'hm3']),
#     expand('data/{dataset}/df_beta_good_lassosum.rds', dataset=[ 'chris',  'hm3']),
#     expand('results/{dataset}/params_grid_lassosum.rds', dataset=['chris', 'hm3'])

# rule all_blocks:
#   input:
#     expand('results/{dataset}/heritability_ldpred2_blocks.rds', dataset=['chris']),
#     expand('results/{dataset}/final_mod_ldpred2_blocks.rds', dataset=['chris']),
#     expand('results/{dataset}/beta_auto_ldpred2_blocks.rds', dataset=['chris']),
#     expand('data/{dataset}/df_beta_good_with_ld_blocks.rds', dataset=['chris']),
#     expand('results/{dataset}/beta_grid_ldpred2_blocks.rds', dataset=['chris'])

# rule all_corr:
#   input:
#     expand('data/ld_ref/cor_{dataset}.rds', dataset=['hm3', 'chris', 'hm3plus'])

# rule corr_hm3plus:
#   input:
#     ldfiles = expand('data/hm3_plus_ref/LD_with_blocks_chr{cc}.rds', cc=range(1,23)),
#     map_rds = 'data/geno/hm3/qc_geno_all_map_gendist.rds'
#   output:
#     corfile = 'data/ld_ref/cor_hm3plus.rds'
#   params:
#     map_hm3 = 'data/map_hm3_plus.rds'
#   threads: 8
#   resources:
#     mem_mb=32000,
#   script:
#     'bin/create_cor.R'

# rule corr:
#   input:
#     ldfiles = expand('data/ld_ref/{{dataset}}/ld_chr{cc}.rds', cc=range(1,23)),
#     map_rds = 'data/geno/{dataset}/qc_geno_all_map_gendist.rds'
#   output:
#     corfile = 'data/ld_ref/cor_{dataset}.rds'
#   params:
#     map_hm3 = ''
#   threads: 8
#   resources:
#     mem_mb=64000,
#     tmpdir='tmp-data'
#   script:
#     'bin/create_cor.R'

rule prs_LDpred2:
  input:
    # train_data = 'data/train_pheno_file.rds',
    map_rds = "data/geno/qc_geno_all_map.rds",
    gwas_rds="results/{pheno}/gwas.rds",
    genotype_rds = "data/geno/qc_geno_all.rds",
    map_ld_rds = ancient(hm3map),
    corfiles = ancient(hm3corr)
    # map_ld_rds = "resources/ld_ref/map_hm3_plus.rds",
    # corfiles = expand("resources/ld_ref/ldref_hm3_plus/LD_with_blocks_chr{chrom}.rds", 
    #                  chrom=range(1,23))
  output:
    "results/{pheno}/ldpred2/heritability_ldpred2.rds",
    "results/{pheno}/ldpred2/final_mod_ldpred2.rds",
    "results/{pheno}/ldpred2/beta_auto_ldpred2.rds",
    "results/{pheno}/ldpred2/df_beta_good_with_ld.rds",
    "results/{pheno}/ldpred2/beta_grid_ldpred2.rds",
    "results/{pheno}/ldpred2/pred_grid_ldpred2.rds",
    "results/{pheno}/ldpred2/params_grid_ldpred2.rds"
  resources:
    mem_mb = 24000,
    tmpdir = "tmp-data"
  threads: 12
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_PRS_ldpred2.R"

rule predict_LDpred2:
  input:
    # train_data = 'data/train_pheno_file.rds',
    # map_rds = 'data/geno/rds/qc_geno_all_map.rds',
    map_rds = "data/geno/qc_geno_all_map.rds",
    gwas_rds="results/{pheno}/gwas.rds",
    map_good = "results/{pheno}/ldpred2/df_beta_good_with_ld.rds",
    genotype_rds = "data/geno/qc_geno_all.rds",
    beta_auto = "results/{pheno}/ldpred2/beta_auto_ldpred2.rds"
  output:
    pred_rds = "results/{pheno}/ldpred2/prs.rds",
    pred_csv = "results/{pheno}/ldpred2/prs.csv",
    map_file = "results/{pheno}/ldpred2/map_prs.rds"
  threads: 16
  resources:
    mem_mb = 32000
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/predict_LDpred2.R"


rule prs_lassosum2:
  message:
    "Run PRS estimation with lassosum2"
  input:
    map_rds = "data/geno/qc_geno_all_map.rds",
    gwas_rds="results/{pheno}/gwas.rds",
    genotype_rds = "data/geno/qc_geno_all.rds",
    map_ld_rds = ancient(hm3map),
    corfiles = ancient(hm3corr)
    # map_ld_rds = "resources/ld_ref/map_hm3_plus.rds", 
    # corfiles = expand("resources/ld_ref/ldref_hm3_plus/LD_with_blocks_chr{chrom}.rds", 
                     # chrom=range(1,23))
  output:
    "results/{pheno}/lassosum2/beta_lassosum2.rds",
    "results/{pheno}/lassosum2/df_beta_good_lassosum.rds",
    "results/{pheno}/lassosum2/params_grid_lassosum.rds"
  threads: 8
  resources:
    mem_mb = 64000,
    tmpdir = "tmp-data"
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_lassosum2.R"


rule predict_lassosum2:
  message:
    "Predict PRS with lassosum2"
  input:
    genotype_rds = "data/geno/qc_geno_all.rds",
    df_beta = "results/{pheno}/lassosum2/df_beta_good_lassosum.rds",
    beta_lassosum = "results/{pheno}/lassosum2/beta_lassosum2.rds",
  output:
    # pred = "results/{pheno}/lassosum2/pred_lassosum2.rds"
    pred_rds = "results/{pheno}/lassosum2/prs.rds",
    pred_csv = "results/{pheno}/lassosum2/prs.csv",
    map_file = "results/{pheno}/lassosum2/map_prs.rds"
  resources:
    mem_mb=14000
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/predict_lassosum2.R"
