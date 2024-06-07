import json

jsonfile = config["gwas_manifest"]
gwases = json.load(open(jsonfile, "r"))
gwas_traits = gwases.keys()

rule create_gwas_rds:
  message: 
    "Creating RDS file for {input}."
  input:
    # gwas_file='data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz',
    gwas_file = lookup("{pheno}/path", within=gwases)
  output:
    gwas_rds = "results/{pheno}/gwas.rds"
  params:
    col_dict = lookup("{pheno}/columns", within=gwases),
    trait_type = lookup("{pheno}/trait_type", within=gwases) 
  resources:
    mem_mb=16000
  script:
    "../scripts/create_gwas_rds.R"
