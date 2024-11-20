import os
from liftover import get_lifter
import pandas as pd

# ld_ref_path = os.path.join


# def main(basepath=".", dataset="ukbb", population="eur"):
def main(input_file="snpinfo_ukbb_hm3", outfile=None):
    
    """Main function to liftover ldblk reference for PRScs method

    Args:
        basepath: string or Path-like object where the ldblk dataset has been downloaded.
        dataset: string, the dataset used for ld reference (default: "ukbb")
        population: string, ethnicity of the reference population (default: "eur")
    Returns:
        a pandas.DataFrame containing the lifted position for the available variants.

        If the variant failed to be converted, then the position is set to "-1" 
        chromosome set to "".
        
    """
    
    # # Create path
    # ldref_dir = f"ldblk_{dataset.lower()}_{population.lower()}"
    # ldref_dir = os.path.join(basepath, ldref_dir)

    # # Start converter
    converter = get_lifter('hg19', 'hg38', one_based=True)

    # # Read snp info
    # snpinfo_file = f"snpinfo_{dataset.lower()}_hm3"
    # snpinfo_file = os.path.join(ldref_dir, snpinfo_file)
    snpinfo = pd.read_csv(input_file, sep="\t")

    # Add column to keep conversion to hg38
    snpinfo["CHR_hg38"] = ""
    snpinfo["BP_hg38"] = -1

    i = 0
    for pos, row in snpinfo.iterrows():
        cc = converter[str(row["CHR"])][row["BP"]]
        if i % 1000 == 0:
            print(f"Computing row: {i}", end="\r")
        if cc:
            snpinfo.loc[pos, "CHR_hg38"] = cc[0][0]
            snpinfo.loc[pos, "BP_hg38"] = cc[0][1]
        i += 1

    if outfile is None:
        outfile = input_file + "_hg38"
    
    print(f"Saving to disk: {outfile}")
    snpinfo.to_csv(outfile, index=False, sep="\t",
                   columns = ["CHR_hg38", 
                              "SNP",
                              "BP_hg38", 
                              "A1", "A2",
                              "MAF"
                              ],
                   header = ["CHR", "SNP", "BP",
                             "A1", "A2", "MAF"])

    return snpinfo


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--basedir", default=".")
    parser.add_argument("-d", "--dataset", default="ukbb", required=False)
    parser.add_argument("-p", "--population", default="eur", required=False)

    if snakemake.input["snp_info"]:
        input_file = snakemake.input["snp_info"]
        output_file = snakemake.output["snp_info_out"]
    else:
        args = parser.parse_args()
        ldref_dir = f"ldblk_{args.dataset.lower()}_{args.population.lower()}"
        input_file = os.path.join(args.basedir, ldref_dir, f"snpinfo_{args.dataset.lower()}_hm3")
        output_file = input_file + "_hg38"

    snpinfo = main(input_file = input_file, outfile = output_file)
    
    print(snpinfo)
    notconv = (snpinfo["BP_hg38"] < 0).sum()
    print(f"Not converted: {notconv} ({notconv/snpinfo.shape[0]}%)")
