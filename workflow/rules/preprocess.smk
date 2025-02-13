import os

wildcard_constraints:
  chrom=r"\d+"

# rule extract_maf_r2:
#   input:
#     vcf = os.path.join(genopath, 'chr{chrom}.dose.vcf.gz')
#   output:
#     snplist = 'data/geno/qc/chr{chrom}_maf_r2.txt'
#   threads: 1
#   resources:
#     mem_mb=32000
#   shell:
#     """
#     zcat {input} | grep -v -e "#" | cut -f 1-8 > {output}
#     """

rule divide_by_chrom:
  output:
    bed = os.path.join(geno_dir, 'geno_chr{chrom}.bed'),
    bim = os.path.join(geno_dir, 'geno_chr{chrom}.bim'),
    fam = os.path.join(geno_dir, 'geno_chr{chrom}.fam')
  conda:
    '../envs/plink.yaml'
  params:
    genotype_conf = genotype_conf
  script:
    '../scripts/divide_by_chrom.py'

rule create_snplist:
  input:
    mafr2 = os.path.join(geno_dir, 'qc', 'chr{chrom}_maf_r2.txt')
  output:
    snplist = os.path.join(geno_dir, 'qc', 'snplist_chr{chrom}.txt')
  resources:
    mem_mb = 8000
  threads: 1
  params:
    maf=0.05,
    r2=0.9
  script:
    'scripts/create_snplist.py'

rule genotype_QC:
  input: 
    unpack(get_reference)
  output:
    bed = os.path.join(geno_dir, 'qc_geno_chr{chrom}.bed'),
    bim = os.path.join(geno_dir, 'qc_geno_chr{chrom}.bim'),
    fam = os.path.join(geno_dir, 'qc_geno_chr{chrom}.fam')
  params:
    prefixi = lambda wildcards, input: input.bed.replace(".bed", ""),
    prefixo = lambda wildcards, output: output.bed.replace(".bed", ""),
    maf = config["maf"],
    hwe = config["hwe"],
    geno = config["genorate"],
    mind = config["mind"]
  resources:
    mem_mb=32000
  conda:
    "../envs/plink.yaml"
  shell:
    """
    plink --bfile {params.prefixi} \
      --memory {resources.mem_mb} \
      --maf {params.maf} \
      --hwe {params.hwe} \
      --geno {params.geno} \
      --mind {params.mind}  \
      --make-bed \
      --out {params.prefixo}
    """

rule write_genotype_mergelist:
  message: "Create plink merge list"
  input:
      # get_all_ref
      get_reference2
  output:
      os.path.join(geno_dir, 'merge_list_all.txt')
  run:
    with open(output[0], 'w') as f:
        for i in input:
            if i.endswith(".bed"):
                fi = i.replace(".bed", "")
                f.write(fi + "\n")

rule merge_all_genotypes:
  input:
    rules.write_genotype_mergelist.output #'merge_list_all.txt'
  output:
    bed = os.path.join(geno_dir, "qc_geno_chrall.bed"),
    bim = os.path.join(geno_dir, "qc_geno_chrall.bim"),
    fam = os.path.join(geno_dir, "qc_geno_chrall.fam")
  params:
    prefixo=lambda wildcards, output: output.bed.replace('.bed', '')
  message:
    "Merge plink files"
  resources:
    mem_mb=72000
  conda:
    "../envs/plink.yaml"
  script:
    "../scripts/plink_merge.py"


rule import_geno_into_r_bychr:
  message:
    "Import genotype datya into R using bigsnpr structure"
  input:
    bed = ancient(os.path.join(geno_dir, "qc_geno_chr{chrom}.bed")),
    bim = ancient(os.path.join(geno_dir, "qc_geno_chr{chrom}.bim")),
    fam = ancient(os.path.join(geno_dir, "qc_geno_chr{chrom}.fam"))
  output:
    os.path.join(geno_dir, "qc_geno_chr{chrom}.rds"),
    os.path.join(geno_dir, "qc_geno_chr{chrom}.bk"),
    os.path.join(geno_dir, "qc_geno_chr{chrom}_map.rds")
  threads: 16 
  resources: 
    mem_mb = 64000,
    tmpdir = "tmp-data"
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/load_genotype_all.R"


rule import_genotype_into_r:
  message:
    "Import genotype datya into R using bigsnpr structure"
  input:
    bed = ancient(os.path.join(geno_dir, "qc_geno_chrall.bed")),
    bim = ancient(os.path.join(geno_dir, "qc_geno_chrall.bim")),
    fam = ancient(os.path.join(geno_dir, "qc_geno_chrall.fam"))
  output:
    os.path.join(geno_dir, "qc_geno_all.rds"),
    os.path.join(geno_dir, "qc_geno_all.bk"),
    os.path.join(geno_dir, "qc_geno_all_map.rds")
  threads: 16 
  resources: 
    mem_mb = 64000,
    tmpdir = "tmp-data"
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/load_genotype_all.R"

rule compute_distance:
  input:
    mapfile = os.path.join(geno_dir, "qc_geno_chr{chrom}_map.rds")
  output:
    os.path.join(geno_dir, "qc_geno_chr{chrom}_map_gendist.rds")
  threads: 12
  resources:
    mem_mb = 32000,
    tmpdir = "resources"
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/compute_genetic_dist.R"


# rule corr_LD_bychr:
#     input:
#         bed = 'data/geno/{dataset}/qc_geno_chr{chrom}.bed',
#         bim = 'data/geno/{dataset}/qc_geno_chr{chrom}.bim',
#         fam = 'data/geno/{dataset}/qc_geno_chr{chrom}.fam',
#         genotype_rds = 'data/geno/{dataset}/qc_geno_all.rds',
#         map_rds = 'data/geno/{dataset}/qc_geno_all_map_gendist.rds'
#     output:
#         cormat = 'data/ld_ref/{dataset}/ld_chr{chrom}.rds',
#         # indepmat = 'ld_ref/ld_indep_blocks_{chrom}.rds'
#     threads: 6,
#     resources:
#         mem_mb = 128000,
#         tmpdir = 'tmp-data'
#     script:
#         'bin/compute_LD_bychr.R'


# rule create_split_onLD:
#   input:
#     ldref = 'data/ld_ref/chris/ld_chr{chrom}.rds'
#   output:
#     ldblocks = 'data/ld_ref_blocks/chris/ld_chr{chrom}.rds',
#     best_split = 'data/ld_ref_blocks/chris/best_split_{chrom}.rds'
#   resources:
#     mem_mb = 64000
#   script:
#     'bin/split_ld_blocks.R'


# rule splits:
#   input:
#     expand('data/ld_ref_blocks/chris/ld_chr{chrom}.rds', chrom=range(1,23))

rule qc_plot:
  input:
    map_rds = os.path.join(geno_dir, "qc_geno_all_map.rds"),
    gwas_rds = os.path.join(data_dir, "gwas", "{pheno}_overall.rds")
  output:
    plot_file = os.path.join(odir, "bad_variants.png"),
    plot_file2 = os.path.join(odir, "beta_distribution.png")
  resources:
    mem_mb=24000
  conda:
    "../envs/bigsnpr.yaml"
  script:
    "../scripts/qc_plot.R"

## rule indep_LD_blocks:
##     input:
##         'ld_ref/ld_chr{chrom}.rds'
##     output:
##         'ld_ref/ld_indep_blocks_{chrom}.rds'
##     resources:
##         mem_mb = 50000,
##         tmpdir = 'tmp-data'
##     script:
##         'bin/indep_LD_blocks.R'

rule create_ldref:
    input:
        ldfiles = expand('ld_ref/ld_chr{cc}.rds', cc=range(1,23)),
        genotype_rds = 'plinkFiles/chrall.rds'
    resources:
        threads = 1,
        tmpdir = "tmp-data"
    output:
        'ld_ref/map_ldref.rds'
    script:
        '../bin/export_ldref.R'

# rule merge_gwas_genotype:
#     input:
#         geno_data = 'plinkFiles/chrall.rds',
#         gwas='data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.rds'
#     resources:
#         threads=12,
#         mem_mb=32000
#     output:
#         dfbeta = 'data/merge_df_beta.rds'
#     script:


rule check_ld_ref:
  input:
    genotype_rds = 'data/geno/{dataset}/qc_geno_all.rds',
    map_rds = 'data/geno/{dataset}/qc_geno_all_map_gendist.rds',
    gwas_rds='data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.rds',
    cormat = 'data/ld_ref/{dataset}/ld_chr{chrom}.rds'
  output:
    out_rds = 'data/ld_ref/{dataset}/condz_in_chr{chrom}.rds',
    out_plot = 'data/ld_ref/{dataset}/condz_in_chr{chrom}.png',
    out_lmb = 'data/ld_ref/{dataset}/condz_in_lambda_chr{chrom}.txt'
  resources:
    mem_mb = 24000
  script:
    "../scripts/plot_LD_by_chr.R"


rule get_rsid:
  message:
    "Download rsid database"
  output:
    rsid_file = protected(rsidfilevar)
  params:
    rsfile = lambda wildcards, output: os.path.basename(output.rsid_file).replace("index.gz", "tsv.gz"),
    resource_dir = resource_dir
  conda:
    "../envs/bcftools.yaml"
  resources:
    mem_mb=8000
  shell:
    """
    if [ ! -f {params.resource_dir}/{params.rsfile} ]
    then
      wget -O {params.resource_dir}/{params.rsfile} https://resources.pheweb.org/{params.rsfile}
    fi
    echo -e "CHROM\tPOS\tRSID\tREF\tALT" | bgzip -c > {output.rsid_file}
    zcat {params.resource_dir}/{params.rsfile} | bgzip -c >> {output.rsid_file}
    tabix -s1 -b2 -e2 -S1 {output.rsid_file}
    """

rule get_ld_ref:
  output:
    hm3plus = protected(hm3map)
  shell:
    """
    wget -O {output.hm3plus} https://figshare.com/ndownloader/files/37802721
    """

rule get_ld_ref_mat:
  message:
    "Download hm3plus correlation matrix"
  output:
    hm3_mat = protected(hm3corr)
  params:
    zipfile = hm3zip,
    outpath = hm3matout
  resources:
    mem_mb = 8000
  shell:
    """
    if [ ! -f {params.zipfile} ]
    then
      gdown -O {params.zipfile} https://drive.google.com/uc?id=17dyKGA2PZjMsivlYb_AjDmuZjM1RIGvs
    fi
    unzip {params.zipfile} -d {params.outpath}
    """

rule get_gwas_formats:
  output:
    protected(get_formatbooks())
  params:
    resource_dir= resource_dir
  shell:
    """
    git clone https://github.com/Cloufield/formatbook.git {params.resource_dir}/formatbook
    """