import json
import os
import itertools as it
from snakemake.utils import validate
from snakemake.rules import expand

# Validate the config file
validate(config, schema="../schemas/config.schema.yaml")

# Genotype configuation
genotype_conf = json.load(open(config["genotype_json"], "r"))
validate(genotype_conf, schema="../schemas/genotype_json.schema.yaml")

# Set the results directory
odir = os.path.join(config["output_dir"], "{pheno}")
data_dir = config["data_dir"]
geno_dir = os.path.join(data_dir, "geno")

# Create directory if they do not exists
# os.makedirs(config["output_dir", exist_ok=True)
# os.makedirs(data_dir, exist_ok=True)
# os.makedirs(geno_dir, exist_ok=True)

# Create resource directory and check existence
resource_dir = config["cache_dir"]
if not os.path.isdir(resource_dir):
  resource_dir = "resources"

# GWAS configuration
jsonfile = config["gwas_manifest"]
gwases = json.load(open(jsonfile, "r"))
gwas_traits = gwases.keys()

# Validate GWAS
# Check build
# TODO: create a configuration file to test different build between genome and GWAS sumstat
for k, v in gwases.items():
  validate(v, schema="../schemas/gwas_manifest.schema.yaml")
  if v["genome_build"] != genotype_conf["build"]:
    raise TypeError("Different genome build between gwas and genotypes. Consider lifting over one of the datasets")


# Create path for rsid file
rsidfilevar = os.path.join(resource_dir, "rsid", f"rsids-v154-{genotype_conf['build']}.index.gz")

# Get hm3 high quality variants from web resources
hm3path = os.path.join(resource_dir, "ld_ref")
hm3zip = os.path.join(hm3path, "ldref_hm3_plus.zip")
hm3corr = expand(os.path.join(hm3path, "ldref_hm3_plus", "LD_with_blocks_chr{chrom}.rds"), chrom=range(1,23))
hm3map = os.path.join(hm3path, "map_hm3_plus.rds")
hm3matout = os.path.join(hm3path, "ldref_hm3_plus")


# Get LD reference configuration
ldref = config["ld_reference"]
lddata = ldref["data"].lower()
ldpop = ldref["population"].lower()


def add_plink_ext(infile):
  myext = ["bed", "bim", "fam"]
  return {k: infile + f".{k}" for k in myext}

def get_reference(wildcards):
  # if genotype_conf["nchrom"] == 1:
  #   ref_files = genotype_conf["plinkfiles"]["1"]
  # else:
  # file_format = genotype_conf["format"]
  # if genotype_conf["divided_by_chrom"]:
  #   if file_format == "plink":
  #     ref_files = genotype_conf["files"][wildcards.chrom]
  #   elif file_format == "vcf":
  #     ref_files = genotype_conf["files"][wildcards.chrom]
  #     ref_files = os.path.basename(ref_files).replace("vcf.gz", "")
  #     ref_files = os.path.join("tmp-data/geno", ref_files)
  # else:
  #   ref_files = genotype_conf["files"][wildcards.chrom]
  ref_files = os.path.join(geno_dir, "geno_chr{chrom}")
  out_dict = add_plink_ext(ref_files)
  return out_dict

def get_reference2(wildcards):
  nchrom = 22
  flist = expand(os.path.join(geno_dir, "qc_geno_chr{chrom}{ext}"), chrom=range(1, nchrom + 1), 
               ext=[".bed", ".bim", ".fam"])
  return flist


# Get files to combine
def get_all_ref(wildcards):
  nchrom = genotype_conf["nchrom"]
  flist = expand(os.path.join(geno_dir, "qc_geno_chr{chrom}{ext}"), chrom=range(1, nchrom + 1), 
                ext=[".bed", ".bim", ".fam"])
  return flist


def get_ldblk_url(wildcards):
  
  # TODO: check for available population 
  # TODO: check for available dataset

  filename = "ldblk_{data}_{population}.tar.gz"
  filename = filename.format(**{
      "data": lddata,
      "population": ldpop
    })

  baseurl = "https://personal.broadinstitute.org/hhuang/public/PRS-CSx/Reference/"
  baseurl += lddata.upper()
  baseurl += f"/{filename}"

  return baseurl

def get_ldblk_zip():
  # lddata = config["ld_data"].lower()
  # ldpop = config["ld_population"].lower()
  
  mydir = os.path.join(resource_dir, "ldblk_{data}_{population}.tar.gz")
  return mydir.format(**{"data":lddata, "population": ldpop})

def get_ldblk_files():
  # lddata = config["ld_data"].lower()
  # ldpop = config["ld_population"].lower()
  # nchroms = genotype_conf["nchrom"]
  
  mydir = os.path.join(resource_dir, "ldblk_{data}_{population}", "ldblk_{data}_chr{chrom}.hdf5")
  return expand(mydir, data=[lddata], population=[ldpop], 
                chrom=range(1, 22))   

def get_ldblk_dir():
  # lddata = config["ld_data"].lower()
  # ldpop = config["ld_population"].lower()
  
  mydir = os.path.join(resource_dir, "ldblk_{data}_{population}")
  return mydir.format(**{"data": lddata, "population": ldpop})

def prscs_beta_collect(wildcards):
  return expand(os.path.join(
    config["output_dir"], "{{pheno}}", 
    "prscs/_pst_eff_a1_b0.5_phiauto_chr{chrom}.txt"), 
    chrom=range(1,genotype_conf["nchrom"] + 1))

# Preprocessing
def target_rule_preproc():
  """Define target rule for preprocessing
  """
  genofiles = [f"{geno_dir}/qc_geno_all.rds", f"{geno_dir}/qc_geno_all.bk",
    f"{geno_dir}/qc_geno_all_map.rds"]
  return genofiles


def expand_chrom(myfile):
  myfile_masked = myfile.replace("pheno", "{pheno}")
  return expand(myfile_masked, chrom=range(1, 23))


def target_rule_preproc_bychr():
  genofiles = expand_chrom(os.path.join(geno_dir, "qc_geno_chr{chrom}.rds"))
  genofiles.append(expand_chrom(os.path.join(geno_dir, "qc_geno_chr{chrom}.bk")))
  genofiles.append(expand_chrom(os.path.join(geno_dir, "qc_geno_chr{chrom}_map.rds")))
  genofiles.append(expand_chrom(os.path.join(geno_dir, "qc_geno_chr{chrom}_map_gendist.rds")))

  return genofiles


# SCT PRS
def target_rule_sct():
  return {"clump": expand(os.path.join(odir, "sct/clump_res_ct.rds"), pheno=gwas_traits),
   "multi_rds": expand(os.path.join(odir, "sct/multi_prs_ct.rds"), pheno=gwas_traits),
   "multi_bk": expand(os.path.join(odir, "sct/multi_prs_ct.bk"), pheno=gwas_traits)}
 
# GWAS
def target_rule_gwas():
  return expand(os.path.join(odir, "gwas.rds"), pheno=gwas_traits)

# Define target rule for computing the PRSs
def target_rule_prs():
  output_files = []
  for k, v in config["prs_algorithms"].items():
    try:
      if v["activate"]:
        alg_path = os.path.join(odir, k.lower())
        preds = expand(os.path.join(alg_path, "prs{ext}"),
                      pheno=gwas_traits, ext=[".rds", ".csv"])
        output_files.extend(preds)
        # if k != "sct":
        #   maps = expand(os.path.join(alg_path, "map_prs.rds"),
        #                 pheno=gwas_traits)
        #   output_files.extend(maps)
    except KeyError:
      pass
  return output_files


def target_rule_plots():
  ofiles = [expand(os.path.join(odir, "bad_variants.png"), pheno=gwas_traits),
  expand(os.path.join(odir, "beta_distribution.png"), pheno=gwas_traits)]
  return ofiles
 

def get_formatbooks():
  return os.path.join(resource_dir, "formatbook", "formatbook.json")

def lift_ld_ref():
  return os.path.join(get_ldblk_dir(), f"snpinfo_{lddata}_hm3_hg38")

# def tmp_ld():
#   alg = "lassosum2"
#   alg_path = os.path.join(odir, alg)
#   ff = expand(os.path.join(alg_path, "beta_auto_chr{chrom}_ldpred2.rds"), 
#   chrom=range(1,23), 
#   pheno=gwas_traits)
#   print(ff)
#   return ff


# # PRS-CS
# def target_rule_prscs():

#   return {"pred_prscs": expand("results/{pheno}/prscs/prs{ext}", 
#                                pheno=gwas_traits, ext=[".rds", ".csv"]),
#           "map_prscs": expand("results/{pheno}/prscs/map_prs.rds", pheno=gwas_traits)}

# # LDPred2
# def target_rule_ldpred2():
#   return {"pred_ldpred2": expand("results/{pheno}/ldpred2/prs{ext}", pheno=gwas_traits,
#                                 ext=[".rds", ".csv"])}


