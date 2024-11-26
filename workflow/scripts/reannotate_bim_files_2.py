import subprocess as sb
import itertools as it
import pandas as pd


def main(bimfile, rsidfile, chrom, outfile="out_annot.bim"):
    """Annotate bim file with RSID.
    Use the strategy of matching chromosome, position and alleles

    """
    bimdf = pd.read_csv(bimfile, delimiter='\t',
                        names=['CHROM', 'ID', 'CM', 'POS', 'A0', 'A1'])

    # Create temporary file to get the tabix out from rsid and from
    # bim file.
    bimtmp = bimfile.replace(".bim", "_reg.bim")
    tabixtmp = bimfile.replace(".bim", "_regout.bim")

    # Extract regions from bim file...
    bimdf.loc[:, ["CHROM", "POS"]].to_csv(bimtmp, index=False, header=False,
                                          sep="\t")
    # Run tabix
    f = open(tabixtmp, "w")
    proc_out = sb.run(["tabix", rsidfile, "-R", bimtmp],
                      stdout=f, text=True)
    f.close()

    if proc_out.returncode == 0:
        regions_rsids = pd.read_csv(tabixtmp, header=0, sep="\t",
                                    names=["CHROM", "POS", "RSID", "A0", "A1"])
    else:
        raise IOError(f"No such file or directory {tabixtmp}!")
    
    # Output file
    bimoutfile = outfile
    bimout = open(bimoutfile, 'w')

    for j, r in enumerate(bimdf.iterrows()):
        pp = r[1]['POS']
        try:
            rs = regions_rsids.loc[pp, :]
            rsarr = rs.to_numpy()
            if len(rsarr.shape) == 1:
                rsarr = rsarr.reshape((1, rsarr.shape[0]))
            rid = compare_alleles(r, rsarr)
        except KeyError:
            rid = r[1]["ID"]

        # Create line to write in the output
        line = [r[1]['CHROM'], rid, r[1]['CM'], r[1]['POS'], r[1]['A0'],
                r[1]['A1']]
        line = [str(l) for l in line]

        # Write to output
        bimout.write('\t'.join(line))
        bimout.write('\n')


def compare_alleles(bim, anno):
    """Function to compare alleles eventually 
    allowing reverse strand and flipping
    """
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
    #
    # main(bimfile=args.bimfile, rsidfile=args.rsidfile, chrom=args.chrom,
    #      outfile=args.out)

    bimfile = snakemake.input.bim
    rsidfile = snakemake.input.rsidfile
    chrom = snakemake.wildcards.chrom
    outfile = snakemake.output.bim_anno
    main(bimfile=bimfile, rsidfile=rsidfile, chrom=chrom, outfile=outfile)
