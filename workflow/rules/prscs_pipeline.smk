# include: "preprocess_forPRS_debug.snakefile"

datageno = 'data/geno/chris'  
ldpath='/scratch/mfilosi/reference/prscs_ld_ref'

rule all_p:
  input:
    pred_file = "results/{pheno}/prscs/pred_prscs.rds",
    map_file = "results/{pheno}/prscs/map_prscs.rds"
    # expand(os.path.join(datageno, 'qc_geno_chr{chrom}_rsid.bim'), chrom=range(1,23))
    # expand("results/{pheno}_pst_eff_a1_b0.5_phi1e-02_chr{chrom}.txt", pheno=pp, chrom=range(1,23))
    # expand("results/prscs/{pheno}_pst_all.txt", pheno=pp)
    # gwas_prscs = 'data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K_prscs.csv'
    # expand("data/gwas/prscs/{pheno}.tsv", pheno=pp)


rule gwas_for_prscs:
  input:
    gwas_rds = "results/{pheno}/gwas.rds"
  output:
    gwas_prscs = "results/{pheno}/gwas_prscs.csv"
  resources:
    mem_mb=16000
  script:
    "../scripts/gwasrds_to_prcs.R"

rule annotate_bim:
  input:
    bim =  "data/geno/qc_geno_chr{chrom}.bim",
    rsidfile = ancient(f"resources/rsid/rsids-v154-{genotype_conf['build']}.index.gz")
  output:
    bim_anno = 'data/geno/qc_geno_chr{chrom}_rsid.bim'),
  resources:
    mem_mb = 8000
  shell:
    """
    python ../scripts/reannotate_bim_files.py {input.bim} {input.rsidfile} --chrom {wildcards.chrom}
    """

rule run_prscs:
  input:
    gwas_prscs = "results/{pheno}/gwas_prscs.csv",
    bim =  "/data/geno/qc_geno_chr{chrom}_rsid.bim",
    ldref = os.path.join(ldpath, "ldblk_ukbb_eur/ldblk_ukbb_chr22.hdf5")
  output:
    beta_prscs = 'results/{pheno}_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt'
  params:
    tmpdir = 'tmp-data',
    prefixld = lambda wildcards, input: os.path.dirname(input.ldref),
    prefixbim = lambda wildcards, input: input.bim.replace('.bim', ''),
    prefixout = lambda wildcards, output: os.path.dirname(output.beta_prscs) + f'/{pheno}'
  resources:
    mem_mb=16000
  threads: 8
  conda:
    "../envs/prscs.yaml"
  shell:
    """
    python PRScs.py \
    --ref_dir={params.prefixld} \
    --bim_prefix={params.prefixbim} \
    --sst_file={input.gwas_prscs} \
    --n_gwas=500000 \
    --chrom={wildcards.chrom} \
    --out_dir={params.prefixout}
    """

rule move_and_collect:
  input:
    beta_prscs = expand('results/{pheno}_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt',
                        chrom=range(1,genotype_conf["nchrom"] + 1))
  output:
    beta_shrinked = 'results/{pheno}/prscs/beta_all.txt'  
  resources:
    mem_mb = 8000
  shell:
    """
    cat {input.beta_prscs} > {output.beta_shrinked}
    """

rule predict_prscs:
  input:
    bim = expand('/data/geno/qc_geno_chr{chrom}_rsid.bim'), 
                 chrom=range(1,genotype_conf["nchrom"] + 1)),
    beta_shrinked = 'results/{pheno}/prscs/beta_all.txt',
    genotype_rds = 'data/geno/chris/qc_geno_all.rds'
  output:
    pred_file = "results/{pheno}/prscs/pred_prscs.rds",
    map_file = "results/{pheno}/prscs/map_prscs.rds"
  resources:
    mem_mb=12000
  script:
    "bin/predict_prscs.R"


def get_ldblk_files(wildcards):
  lddata = config["ld_data"].lower()
  ldpop = config["ld_population"].lower()
  nchroms = genotype_conf["nchrom"]
  
  mydir = "resources/ldblk_{data}_{population}/ldblk_{data}_{population}_chr{chrom}.hdf5"
  return expand(mydir, data=[lddata], population=[ldpop], 
                chrom=range(1, nchroms + 1))   

def get_ldblk_dir(wildcards):
  lddata = config["ld_data"].lower()
  ldpop = config["ld_population"].lower()
  nchroms = genotype_conf["nchrom"]
  
  mydir = "resources/ldblk_{data}_{population}"
  return mydir.format(**{"data": lddata}, "population": ldpop})

rule get_ld_ref:
  output: 
    dir(get_ldblk_dir())
  params:
    url = get_ldblk_url()
  shell:
    """
    wget -O {output} {params.url}
    tar -zxvf {output}
    """

