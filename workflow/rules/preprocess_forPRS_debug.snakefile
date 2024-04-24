genopath="/scratch/compgen/data/genetics/CHRIS13K/Imputation/HRC"

wildcard_constraints:
    chrom="\d+"

rule all:
  input:
    expand('data/ld_ref/hm3/ld_chr{chrom}.rds', chrom=range(1,23)),

rule alltmp:
  input:
    # expand('data/geno/qc/snplist_chr{cc}.txt', cc=range(1,23))
    expand('data/geno/{dd}/qc_geno_all.rds', dd=['chris', 'hm3']),
    expand('data/geno/{dd}/qc_geno_all.bk', dd=['chris', 'hm3']),
    expand('data/geno/{dd}/qc_geno_all_map.rds', dd=['chris', 'hm3']),
    expand('data/ld_ref/{dataset}/ld_chr{chrom}.rds', dataset=['chris', 'hm3'], chrom=range(1,23)),
    expand('data/{dataset}/bad_variants.png', dataset=['chris', 'hm3']),
    expand('data/{dataset}/beta_distribution.png', dataset=['chris', 'hm3'])
    # 'ld_ref/ld_chr11.rds'
    # 'data/geno/rds/qc_geno_all.rds',
    # 'data/geno/rds/qc_geno_all.bk',
    # 'data/geno/rds/qc_geno_all_map.rds'

rule create_gwas_rds:
  input:
    gwas_file='data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz',
  output:
    # gwas_rds='data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.rds'
    gwas_rds='data/gwas/eGFR_overall.rds'
  message: "Creating RDS file from: {input}."
  resources:
    mem_mb=16000
  script:
    'bin/create_gwas_rds.R'

rule filter_gwas_rds:
  input:
    gwas_rds='data/eGFR_overall_filt.rds'
  output:
    gwas_rds_filt = 'data/gwas/eGFR_overall_filt.rds'
  resources:
    mem_mb=16000
  script:
    'bin/select_loci.R'

# Uncomment the following lines and add to the shell script below if you want to use bcftools
# bcftools is very slow compare to bash scripting as did here.
# FMT_STRING = r"""%CHROM\t%POS\t%REF\t%ALT{0}\t%INFO/R2\t%INFO/MAF\n"""
# bcftools query -f {FMT_STRING:q} -H -o {output} {input}

rule extract_maf_r2:
  input:
    vcf = os.path.join(genopath, 'chr{chrom}.dose.vcf.gz')
  output:
    snplist = 'data/geno/qc/chr{chrom}_maf_r2.txt'
  threads: 1
  resources:
    mem_mb=32000
  shell:
    """
    zcat {input} | grep -v -e "#" | cut -f 1-8 > {output}
    """

rule create_snplist:
  input:
    mafr2 = 'data/geno/qc/chr{chrom}_maf_r2.txt'
  output:
    snplist = 'data/geno/qc/snplist_chr{chrom}.txt'
  resources:
    mem_mb = 8000
  threads: 1
  params:
    maf=0.05,
    r2=0.9
  script:
    'bin/create_snplist.py'


rule genotype_QC:
    input: 
        bed=os.path.join(genopath, 'plinkFormat', 'chr{chrom}.bed'), 
        bim=os.path.join(genopath, 'plinkFormat', 'chr{chrom}.bim'), 
        fam=os.path.join(genopath, 'plinkFormat', 'chr{chrom}.fam'),
        snplist = 'data/geno/qc/snplist_chr{chrom}.txt'
    output:
        bed='data/geno/chris/qc_geno_chr{chrom}.bed',
        bim='data/geno/chris/qc_geno_chr{chrom}.bim',
        fam='data/geno/chris/qc_geno_chr{chrom}.fam'
    resources:
        mem_mb=32000
    params:
        prefixi=lambda wildcars, input: input.bed.replace(".bed", ""),
        prefixo=lambda wildcars, output: output.bed.replace(".bed", "")
    shell:
        """
        plink --bfile {params.prefixi} --memory {resources.mem_mb} --maf 0.01 --hwe 1e-50 --geno 0.1 --mind 0.1 --make-bed --out {params.prefixo} --extract {input.snplist}
        """

rule write_genotype_mergelist2:
    message: "Create plink merge list"
    input:
        expand('data/geno/{{dataset}}/qc_geno_chr{cc}{ext}', cc=range(1,23), ext=['.bed', '.bim', '.fam'])
    output:
        'merge_list_{dataset}_all.txt'
    run:
        with open(output[0], 'w') as f:
            for i in input:
                if i.endswith(".bed"):
                    fi = i.replace(".bed", "")
                    f.write(fi + "\n")

rule merge_all_genotypes:
  input:
    'merge_list_{dataset}_all.txt'
  output:
    bed = 'data/geno/{dataset}/qc_geno_chrall.bed',
    bim = 'data/geno/{dataset}/qc_geno_chrall.bim',
    fam = 'data/geno/{dataset}/qc_geno_chrall.fam'
  message: "Merge plink files"
  params:
    prefixo=lambda wildcards, output: output.bed.replace('.bed', '')
  resources:
    mem_mb=32000
  shell:
    'plink --merge-list {input} --memory {resources.mem_mb} --make-bed --out {params.prefixo}'

rule reference_data_genotype:
  input:
    bed = 'data/geno/{dataset}/qc_geno_chrall.bed',
    bim = 'data/geno/{dataset}/qc_geno_chrall.bim',
    fam = 'data/geno/{dataset}/qc_geno_chrall.fam'
  output:
    'data/geno/{dataset}/qc_geno_all.rds',
    'data/geno/{dataset}/qc_geno_all.bk',
    'data/geno/{dataset}/qc_geno_all_map.rds'
  threads: 16
  resources: 
    mem_mb = 64000,
    tmpdir = 'tmp-data'
  script:
    'bin/load_genotype_all.R'

rule merge_hm3_plus:
  input:
    bim = 'data/geno/chris/qc_geno_chr{chrom}.bim',
    map_hm3_plus = "data/map_hm3_plus.rds"
  output:
    snplist_file = 'data/geno/hm3_chris_snplist_chr{chrom}.txt'
  resources:
    mem_mb = 24000
  threads: 1
  script:
    'bin/merge_with_hm3_plus.R'

rule select_snps:
  input:
    bed = 'data/geno/chris/qc_geno_chr{chrom}.bed',
    bim = 'data/geno/chris/qc_geno_chr{chrom}.bim',
    fam = 'data/geno/chris/qc_geno_chr{chrom}.fam',
    snplist_file = 'data/geno/hm3_chris_snplist_chr{chrom}.txt'
  output:
    bed = 'data/geno/hm3/qc_geno_chr{chrom}.bed',
    bim = 'data/geno/hm3/qc_geno_chr{chrom}.bim',
    fam = 'data/geno/hm3/qc_geno_chr{chrom}.fam'
  params:
    prefixi=lambda wildcars, input: input.bed.replace(".bed", ""),
    prefixo=lambda wildcars, output: output.bed.replace(".bed", "")
  resources:
    mem_mb = 12000
  shell:
    'plink --bfile {params.prefixi} --memory {resources.mem_mb} --extract {input.snplist_file} --make-bed --out {params.prefixo}'


rule compute_distance:
  input:
    mapfile = 'data/geno/{dataset}/qc_geno_all_map.rds'
  output:
    'data/geno/{dataset}/qc_geno_all_map_gendist.rds'
  threads: 12
  resources:
    mem_mb = 32000,
    tmpdir = 'tmp-data'
  script:
    'bin/compute_genetic_dist.R'

rule create_pheno_for_chris:
    input:
        phenofile = "data/PRS_eGFR_adj_corrsep.txt",
        fam = 'data/geno/chris/qc_geno_chrall.fam'
        # fam=os.path.join(genopath, 'plinkFormat', 'chr1.fam')
    output:
        'train_sample_for_geno.txt',
        'data/train_pheno_file.rds',
        'data/test_pheno_file.rds',
    script:
        'bin/load_pheno_and_sample_list.R'

rule corr_LD_bychr:
    input:
        bed = 'data/geno/{dataset}/qc_geno_chr{chrom}.bed',
        bim = 'data/geno/{dataset}/qc_geno_chr{chrom}.bim',
        fam = 'data/geno/{dataset}/qc_geno_chr{chrom}.fam',
        genotype_rds = 'data/geno/{dataset}/qc_geno_all.rds',
        map_rds = 'data/geno/{dataset}/qc_geno_all_map_gendist.rds'
    output:
        cormat = 'data/ld_ref/{dataset}/ld_chr{chrom}.rds',
        # indepmat = 'ld_ref/ld_indep_blocks_{chrom}.rds'
    threads: 6,
    resources:
        mem_mb = 128000,
        tmpdir = 'tmp-data'
    script:
        'bin/compute_LD_bychr.R'


rule create_split_onLD:
  input:
    ldref = 'data/ld_ref/chris/ld_chr{chrom}.rds'
  output:
    ldblocks = 'data/ld_ref_blocks/chris/ld_chr{chrom}.rds',
    best_split = 'data/ld_ref_blocks/chris/best_split_{chrom}.rds'
  resources:
    mem_mb = 64000
  script:
    'bin/split_ld_blocks.R'


rule splits:
  input:
    expand('data/ld_ref_blocks/chris/ld_chr{chrom}.rds', chrom=range(1,23))

rule qc_plot:
  input:
    map_rds = 'data/geno/{dataset}/qc_geno_all_map_gendist.rds',
    gwas_rds = 'data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.rds',
  output:
    plot_file = 'data/{dataset}/bad_variants.png',
    plot_file2 = 'data/{dataset}/beta_distribution.png'
  resources:
    mem_mb=24000
  script:
    'bin/qc_plot.R'
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
        'bin/export_ldref.R'


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


rule plot_ld_all:
  input:
    expand('data/ld_ref/{dataset}/condz_in_chr{chrom}.png', dataset="hm3", chrom=range(1,23))

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
    "bin/plot_LD_by_chr.R"

