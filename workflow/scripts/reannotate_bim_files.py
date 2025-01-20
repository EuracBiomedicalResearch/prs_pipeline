import subprocess as sb
import pandas as pd
import itertools as it


def main(bimfile, rsidfile, chrom, outfile="out_annot.bim"):
    # get rsid
    bimdf = pd.read_csv(bimfile, delimiter='\t',
                        names=['CHROM', 'ID', 'CM', 'POS', 'A0', 'A1'])
    # bimoutfile = bimfile.replace(".bim", "_rsid.bim")
    bimoutfile = outfile
    bimout = open(bimoutfile, 'w')
    rsid = ['' for i in range(bimdf.shape[0])]

    for j, r in enumerate(bimdf.iterrows()):
        pp = r[1]['POS']
        cc = sb.run(['tabix', rsidfile, f'{chrom}:{pp}-{pp}'],
                    capture_output=True)
        if cc.stdout:
            rs = cc.stdout.decode()
            rsidslist = [l for l in rs.split('\n') if l != '']
            rsidslist = [l.split('\t') for l in rsidslist]
            if len(rsidslist) > 1:
                # print(r)
                # print(rsidslist)
                rid = compare_alleles(r, rsidslist)
            else:
                rid = rsidslist[0][2]
            if rid:
                rsid[j] = rid
            else:
                rsid[j] = r[1]['ID']
        else:
            rsid[j] = r[1]['ID']

        line = [r[1]['CHROM'], rsid[j], r[1]['CM'], r[1]['POS'], r[1]['A0'],
                r[1]['A1']]
        line = [str(l) for l in line]
        bimout.write('\t'.join(line))
        bimout.write('\n')

    # bimdf['RSID'] = rsid
    return


def compare_alleles(bim, anno):
    bima0 = bim[1]['A0']
    bima1 = bim[1]['A1']
    for a in anno:
        a0 = a[3].split(',')
        a1 = a[4].split(',')
        for g1, g2 in it.product(a0, a1):
            if (bima0 == g1 and bima1 == g2) or (bima1 == g1 and bima0 == g2):
                return a[2]


if __name__ == '__main__':
    # import argparse
    #
    # parser = argparse.ArgumentParser()
    # parser.add_argument("bimfile")
    # parser.add_argument("rsidfile")
    # parser.add_argument("--chrom", type=int, default=22)
    # parser.add_argument("--out", default="out_annot.bim")
    #
    # args = parser.parse_args()
    # main(bimfile=args.bimfile, rsidfile=args.rsidfile, chrom=args.chrom,
    #      outfile=args.out)

    bimfile = snakemake.input.bim
    rsidfile = snakemake.input.rsidfile
    chrom = snakemake.wildcards.chrom
    outfile = snakemake.output.bim_anno
    main(bimfile=bimfile, rsidfile=rsidfile, chrom=chrom, outfile=outfile)
