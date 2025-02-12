import subprocess as sb
import os
import shutil
import json


def main(genotype_conf, chrom, ofile):
    """This function will separate the genetic data into chromosomes
    """
    if genotype_conf["divided_by_chrom"]:
        for k, v in ofile.items():
            infile = genotype_conf["files"][chrom] + f".{k}"
            shutil.copy(infile, v)
    else:
        print(genotype_conf["files"]["1"])
        ofile = ofile["bed"].replace(".bed", "")
        sb.run(["plink",  "--bfile",
                genotype_conf["files"]["1"],
                "--chr", str(chrom),
                "--make-bed",
                "--out", ofile]
               )


if __name__ == "__main__":
    genotype_conf = snakemake.params.genotype_conf
    print(snakemake.output)
    ofile = snakemake.output
    chrom = snakemake.wildcards.chrom
    # ofile = ofile.replace(".bed", "")
    main(genotype_conf=genotype_conf, chrom=chrom, ofile=ofile)
