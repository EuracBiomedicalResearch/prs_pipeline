import json

jsonfile = config["gwas_manifest"]
gwases = json.load(open(jsonfile, "r"))
gwas_traits = gwases.keys()

rule create_gwas_rds:
  message: 
    "Creating RDS file for {input}."
  input:
    # gwas_file='data/19-Dec-22_MW_eGFR_overall_EA_subtracted_CHRIS5K.txt.gz',
    gwas_file = lookup("{pheno}/path", within=gwases),
    formats = get_formatbooks()
  output:
    gwas_rds = os.path.join(odir, "gwas.rds")
  params:
    gwas_conf = lookup("{pheno}", within=gwases),
    res_dir = resource_dir
  resources:
    mem_mb=16000
  script:
    "../scripts/create_gwas_rds.R"
