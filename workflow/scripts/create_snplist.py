import csv
import os

mafthr = float(snakemake.params['maf'])
r2thr = float(snakemake.params['r2'])
# mafthr = 0.05
# r2thr = 0.7

outfile = snakemake.output['snplist']
infile = snakemake.input['mafr2']
# infile = '/scratch/mfilosi/CKDgenPRS/data/geno/qc/chr22_maf_r2.txt'
# outfile = '/scratch/mfilosi/CKDgenPRS/data/geno/qc/snplist_chr22.txt'


def parse_fields(x):
    tmplist = [l.split("=")  for l in x.split(";")]
    mydict = {l[0]:float(l[1]) for l in tmplist if len(l)>1}
    return(mydict)

with open(outfile, 'w') as fout:
    with open(infile, 'r') as fin:
        # header = next(fin)
        for l in fin:
            line = l.strip().split('\t')
            fdict = parse_fields(line[7])
            try:
                if fdict['MAF'] >= mafthr and fdict['R2'] >= r2thr:
                    fout.write(line[2])
                    fout.write('\n')
            except KeyError:
                pass

