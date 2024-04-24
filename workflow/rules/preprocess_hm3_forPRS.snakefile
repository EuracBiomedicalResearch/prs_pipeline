include: "preprocess_forPRS_debug.snakefile"

wildcard_constraints:
    chrom="\d+"

rule allhm3:
  input: 
    expand('data/ld_ref/hm3/ld_chr{cc}.rds', cc=range(1,23))

rule merge_hm3_plus:
  input:
    bim = 'data/geno/plinkFormat/qc_geno_chrall.bim',
    map_hm3_plus = "data/map_hm3_plus.rds"
  output:
    snplist_file = 'data/geno/hm3_chris_snplist.txt'
  resources:
    mem_mb = 24000
  threads: 1
  script:
    'bin/merge_with_hm3_plus.R'

rule select_snps:
  input:
    bed = 'data/geno/chris/qc_geno_chrall.bed',
    bim = 'data/geno/chris/qc_geno_chrall.bim',
    fam = 'data/geno/chris/qc_geno_chrall.fam',
    snplist_file = 'data/geno/hm3_chris_snplist.txt'
  output:
    bed = 'data/geno/hm3/qc_geno_chrall.bed',
    bim = 'data/geno/hm3/qc_geno_chrall.bim',
    fam = 'data/geno/hm3/qc_geno_chrall.fam'
  params:
    prefixi=lambda wildcars, input: input.bed.replace(".bed", ""),
    prefixo=lambda wildcars, output: output.bed.replace(".bed", "")
  resources:
    mem_mb = 32000
  shell:
    'plink --bfile {params.prefixi} --memory {resources.mem_mb} --extract {input.snplist_file} --make-bed --out {params.prefixo}'

# rule plink_to_rds:
#   input:
#     bed = 'data/geno/hm3/qc_geno_chrall.bed',
#     bim = 'data/geno/hm3/qc_geno_chrall.bim',
#     fam = 'data/geno/hm3/qc_geno_chrall.fam'
#   output:
#     'data/geno/hm3/qc_geno_all.rds',
#     'data/geno/hm3/qc_geno_all.bk',
#     'data/geno/hm3/qc_geno_all_map.rds'
#   threads: 8
#   resources: 
#     mem_mb = 64000,
#     tmpdir = 'tmp-data'
#   script:
#     'bin/load_genotype_all.R'

# rule compute_distance_2:
#   input:
#     mapfile = 'data/geno/hm3/qc_geno_all_map.rds'
#   output:
#     'data/geno/hm3/qc_geno_all_map_gendist.rds'
#   threads: 12
#   resources:
#     mem_mb = 32000,
#     tmpdir = 'tmp-data'
#   script:
#     'bin/compute_genetic_dist.R'

# rule corr_LD_bychr2:
#     input:
#         bed = 'data/geno/plinkFormat/qc_geno_chr{chrom}.bed',
#         bim = 'data/geno/plinkFormat/qc_geno_chr{chrom}.bim',
#         fam = 'data/geno/plinkFormat/qc_geno_chr{chrom}.fam',
#         genotype_rds = 'data/geno/hm3/qc_geno_all.rds',
#         map_rds = 'data/geno/hm3/qc_geno_all_map_gendist.rds'
#     output:
#         cormat = 'data/ld_ref/hm3/ld_chr{chrom}.rds',
#         # indepmat = 'ld_ref/ld_indep_blocks_{chrom}.rds'
#     threads: 12,
#     resources:
#         mem_mb = 64000,
#         tmpdir = 'tmp-data'
#     script:
#         'bin/compute_LD_bychr.R'

# rule reference_data_hm3_plus:
#   input:
#     bed = 'data/geno/hm3/geno_hm3_chrall.bed',
#     bim = 'data/geno/hm3/geno_hm3_chrall.bim',
#     fam = 'data/geno/hm3/geno_hm3_chrall.fam'
#   output:
#     'data/geno/hm3/geno_hm3_all.rds',
#     'data/geno/hm3/geno_hm3_all.bk',
#     'data/geno/hm3/geno_hm3_all_map.rds'
#   threads: 8
#   resources: 
#     mem_mb = 64000,
#     tmpdir = 'tmp-data'
#   script:
#     'bin/load_genotype_all.R'

