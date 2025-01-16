import subprocess as sb
import os
import shutil
import json


def main(genotype_conf, chrom, ofile):
    """This function will separate the genetic data into chromosomes
    """
    if genotype_conf["divided_by_chrom"]:
        shutil.copy(genotype_conf["files"][chrom], ofile)
    else:
        print(genotype_conf["files"]["1"])
        sb.run(["plink",  "--bfile",
                genotype_conf["files"]["1"],
                "--chr", str(chrom),
                "--make-bed",
                "--out", ofile]
               )


if __name__ == "__main__":
    genotype_conf = snakemake.params.genotype_conf
    ofile = snakemake.output.bed
    chrom = snakemake.wildcards.chrom
    ofile = ofile.replace(".bed", "")
    main(genotype_conf=genotype_conf, chrom=chrom, ofile=ofile)
