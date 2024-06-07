import json
import os
import itertools as it

genotype_conf = json.load(open(config["genotype_json"], "r"))
jsonfile = config["gwas_manifest"]
gwases = json.load(open(jsonfile, "r"))
gwas_traits = gwases.keys()

# Create resource directory and check existence
resource_dir = config["cache_dir"]
if not os.path.isdir(resource_dir):
  resource_dir = "resources"

# Create path for rsid file
rsidfilevar = os.path.join(resource_dir, "rsid", 
                           f"rsids-v154-{genotype_conf['build']}.index.gz")
# rsidfilevar = f"{resource_dir}/rsid/rsids-v154-{genotype_conf['build']}.index.gz" 


hm3corr = expand(os.path.join(resource_dir, "ld_ref", "ldref_hm3_plus", 
                                    "LD_with_blocks_chr{chrom}.rds"), 
                     chrom=range(1,23))
hm3map = os.path.join(resource_dir, "ld_ref", "map_hm3_plus.rds")


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


def get_ldblk_url(wildcards):
  
  # TODO: check for available population 
  # TODO: check for available dataset

  filename = "ldblk_{data}_{population}.tar.gz"
  filename = filename.format(**{
      "data": config["ld_data"].lower(),
      "population": config["ld_population"].lower()
    })

  baseurl = "https://personal.broadinstitute.org/hhuang/public/PRS-CSx/Reference/"
  baseurl += config["ld_data"].upper()
  baseurl += f"/{filename}"

  return baseurl

def get_ldblk_zip():
  lddata = config["ld_data"].lower()
  ldpop = config["ld_population"].lower()
  nchroms = genotype_conf["nchrom"]
  
  mydir = os.path.join(resource_dir, "ldblk_{data}_{population}.tar.gz")
  return mydir.format(**{"data":lddata, "population": ldpop})

def get_ldblk_files():
  lddata = config["ld_data"].lower()
  ldpop = config["ld_population"].lower()
  nchroms = genotype_conf["nchrom"]
  
  mydir = os.path.join(resource_dir, "ldblk_{data}_{population}", "ldblk_{data}_chr{chrom}.hdf5")
  return expand(mydir, data=[lddata], population=[ldpop], 
                chrom=range(1, nchroms + 1))   

def get_ldblk_dir():
  lddata = config["ld_data"].lower()
  ldpop = config["ld_population"].lower()
  
  mydir = os.path.join(resource_dir, "ldblk_{data}_{population}")
  return mydir.format(**{"data": lddata, "population": ldpop})


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
  return {"clump": expand('results/{pheno}/sct/clump_res_ct.rds', pheno=gwas_traits),
   "multi_rds": expand('results/{pheno}/sct/multi_prs_ct.rds', pheno=gwas_traits),
   "multi_bk": expand('results/{pheno}/sct/multi_prs_ct.bk', pheno=gwas_traits)}
 
# GWAS
def target_rule_gwas():
  return expand("results/{pheno}/gwas.rds", pheno=gwas_traits)

# Define target rule for computing the PRSs
def target_rule_prs():
  out_path = os.path.join("results", "{pheno}")
  output_files = []
  for k, v in config["prs_algorithms"].items():
    try:
      # print(v["activate"])
      if v["activate"]:
        alg_path = os.path.join(out_path, k.lower())
        preds = expand(os.path.join(alg_path, "prs{ext}"),
                      pheno=gwas_traits, ext=[".rds", ".csv"])
        output_files.extend(preds)
        if k != "sct":
          maps = expand(os.path.join(alg_path, "map_prs.rds"),
                        pheno=gwas_traits)
          output_files.extend(maps)
    except KeyError:
      pass
  return output_files

# PRS-CS
def target_rule_prscs():
  # return expand("results/{pheno}/gwas_prscs.csv", pheno=gwas_traits)
  # return expand("results/{pheno}/prscs/beta_all.txt", pheno=gwas_traits)

  return {"pred_prscs": expand("results/{pheno}/prscs/prs{ext}", 
                               pheno=gwas_traits, ext=[".rds", ".csv"]),
          "map_prscs": expand("results/{pheno}/prscs/map_prs.rds", pheno=gwas_traits)}

# LDPred2
def target_rule_ldpred2():
  return {"pred_ldpred2": expand("results/{pheno}/ldpred2/prs{ext}", pheno=gwas_traits,
                                ext=[".rds", ".csv"])}


