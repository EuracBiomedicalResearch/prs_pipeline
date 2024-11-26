# rule all_p:
#   input:
#     pred_file = "results/{pheno}/prscs/pred_prscs.rds",
#     map_file = "results/{pheno}/prscs/map_prscs.rds"

rule gwas_for_prscs:
  input:
    gwas_rds = os.path.join(odir, "gwas.rds")
  output:
    gwas_prscs = os.path.join(odir, "gwas_prscs.csv")
  resources:
    mem_mb = 16000
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/gwasrds_to_prcs.R"

rule annotate_bim:
  input:
    bim =  os.path.join(geno_dir, "qc_geno_chr{chrom}.bim"),
    rsidfile = ancient(rsidfilevar)
  output:
    bim_anno = os.path.join(geno_dir, "qc_geno_chr{chrom}_rsid.bim"),
  conda:
    "../envs/bcftools.yaml"
  resources:
    mem_mb = 8000
  # shell:
  #   """
  #   python scripts/reannotate_bim_files.py {input.bim} {input.rsidfile} --chrom {wildcards.chrom} --out {output.bim_anno}
  #   """
  script:
    "../scripts/reannotate_bim_files_2.py"

rule run_prscs:
  input:
    gwas_prscs = os.path.join(odir, "gwas_prscs.csv"),
    bim =  os.path.join(geno_dir, "qc_geno_chr{chrom}_rsid.bim"),
    ldref = ancient(get_ldblk_files())
  output:
    beta_prscs = os.path.join(odir, "prscs/_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt")
  params:
    tmpdir = "tmp-data",
    # prefixld = lambda wildcards, input: os.path.dirname(input.ldref),
    prefixld = get_ldblk_dir(),
    prefixbim = lambda wildcards, input: input.bim.replace(".bim", ""),
    prefixout = lambda wildcards, output: os.path.dirname(output.beta_prscs) + "/"
  resources:
    mem_mb=16000
  threads: 8
  conda:
    "../envs/prscs.yaml"
  shell:
    """
    PRScs \
    --ref_dir={params.prefixld} \
    --bim_prefix={params.prefixbim} \
    --sst_file={input.gwas_prscs} \
    --n_gwas=500000 \
    --chrom={wildcards.chrom} \
    --out_dir={params.prefixout}
    """

rule move_and_collect:
  input:
    beta_prscs = prscs_beta_collect 
    # beta_prscs = expand("results/{{pheno}}/prscs/_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt",
    #                     chrom=range(1,genotype_conf["nchrom"] + 1))
  output:
    beta_shrinked = os.path.join(odir, "prscs/beta_all.txt")
  resources:
    mem_mb = 8000
  shell:
    """
    cat {input.beta_prscs} > {output.beta_shrinked}
    """

rule predict_prscs:
  input:
    bim = expand(os.path.join(geno_dir, "qc_geno_chr{chrom}_rsid.bim"), 
                 chrom=range(1, genotype_conf["nchrom"] + 1)),
    beta_shrinked = os.path.join(odir, "prscs/beta_all.txt"),
    genotype_rds = os.path.join(geno_dir, "qc_geno_all.rds"),
  output:
    pred_file = os.path.join(odir, "prscs/prs.rds"),
    pred_csv = os.path.join(odir, "prscs/prs.csv"),
    map_file = os.path.join(odir, "prscs/map_prs.rds"),
  params: 
    gwas_conf = lookup("{pheno}", within=gwases),
    ld_dir = get_ldblk_dir(),
    lddata = lddata
  resources:
    mem_mb=12000
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/predict_prscs.R"


rule get_ld_ref_prscs:
  message:
    "Download ld reference"
  output: 
    protected(get_ldblk_files())
  params:
    url = get_ldblk_url,
    odir = get_ldblk_zip()
  resources:
    mem_mb=12000
  shell:
    """
    if [ ! -f {params.odir} ];
    then
      wget -O {params.odir} {params.url}
    fi
    tar -zxvf {params.odir} -C resources
    """

rule lift_reference:
  message:
    "Lift LD reference for PRScs"
  input:
    snp_info = os.path.join(get_ldblk_dir(), f"snpinfo{{lddata}}_hm3")
  output:
    snp_info_out = os.path.join(get_ldblk_dir(), f"snpinfo{{lddata}}_hm3_hg38")
  conda:
    "../envs/prscs.yaml"
  script:
    "../scripts/lift_ld_ref.py"
