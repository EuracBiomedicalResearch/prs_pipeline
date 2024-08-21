import subprocess as sb
import shutil

merge_list = []
with open(snakemake.input[0], "r") as f:
    for ll in f.readlines():
        merge_list.append(ll.strip())
nfiles = len(merge_list)

if nfiles == 1:
    for ext in ["bed", "bim", "fam"]:
        shutil.copy2(f"{merge_list[0]}.{ext}", snakemake.output[ext])
else:
    sb.run(["plink", "--merge-list", snakemake.input[0],
            "--memory", str(snakemake.resources["mem_mb"]) ,
            "--make-bed", "--out", snakemake.params["prefixo"]])
