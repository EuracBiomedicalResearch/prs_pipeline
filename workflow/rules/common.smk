import json
import itertools as it

genotype_conf = json.load(open(config["genotype_json"], "r"))
jsonfile = config["gwas_manifest"]
gwases = json.load(open(jsonfile, "r"))
gwas_traits = gwases.keys()

def add_plink_ext(infile):
  myext = ["bed", "bim", "fam"]
  return {k: infile + f".{k}" for k in myext}

def get_reference(wildcards):
  if genotype_conf["nchrom"] == 1:
    ref_files = genotype_conf["plinkfiles"]["1"]
  else:
    ref_files = genotype_conf["plinkfiles"][wildcards.chrom]

  out_dict = add_plink_ext(ref_files)
  return out_dict

# Get files to combine
def get_all_ref(wildcards):
  flist = []
  for i, e in it.product(range(genotype_conf["nchrom"]), 
                         ['.bed', '.bim', '.fam']):
    chrom = i + 1
    flist.append('data/geno/qc_geno_chr{chrom}{e}')

  nchrom = genotype_conf["nchrom"]
  flist = expand('data/geno/qc_geno_chr{chrom}{ext}', chrom=range(1, nchrom + 1), 
                ext=['.bed', '.bim', '.fam'])
  return flist

# Preprocessing
def target_rule_preproc():
  """Define target rule for preprocessing
  """
  genofiles = ['data/geno/qc_geno_all.rds', 'data/geno/qc_geno_all.bk',
    'data/geno/qc_geno_all_map.rds']
  return genofiles


def expand_chrom(myfile):
  return expand(myfile, chrom=range(1, genotype_conf["nchrom"] + 1))


def target_rule_preproc_bychr():
  genofiles = expand_chrom("data/geno/qc_geno_chr{chrom}.rds")
  genofiles.append(expand_chrom("data/geno/qc_geno_chr{chrom}.bk"))
  genofiles.append(expand_chrom("data/geno/qc_geno_chr{chrom}_map.rds"))
  genofiles.append(expand_chrom("data/geno/qc_geno_chr{chrom}_map_gendist.rds"))
    
  return genofiles


# SCT PRS
def target_rule_sct():
  return {"clump": expand('results/{pheno}/clump_res_ct.rds', pheno=gwas_traits),
   "multi_rds": expand('results/{pheno}/multi_prs_ct.rds', pheno=gwas_traits),
   "multi_bk": expand('results/{pheno}/multi_prs_ct.bk', pheno=gwas_traits)}
 

# GWAS
def target_rule_gwas():
  return expand("results/{pheno}/gwas.rds", pheno=gwas_traits)


def get_ldblk_url():
  
  # TODO: check for available population 
  # TODO: check for available dataset

  filename = "ldblk_{data}_{population}.tar.gz"
  filename.format(**{
      "data": config["ld_data"].upper(),
      "population": config["ld_population"].lower()
    })


  baseurl = "https://personal.broadinstitute.org/hhuang/public/PRS-CSx/Reference/"
  baseurl += config["ld_data"].upper()
  baseurl += f"/{filename}"

  return baseurl

# def get_ldblk_dir():
#   mydir = "resources/ld{}"