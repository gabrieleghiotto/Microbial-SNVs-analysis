import os
import argparse
import pandas as pd

def load_snv_files(input_dir):
    snv = pd.DataFrame()
    for subdir, _, files in os.walk(input_dir):
        for file in files:
            if 'SNVs' in file:
                labels = file.replace('.IS_SNVs.tsv', '')
                df = pd.read_csv(os.path.join(subdir, file), sep="\t")
                df["experiment"] = labels
                snv = pd.concat([snv, df], ignore_index=True)
    return snv

def load_scaffold_info(input_dir):
    scaf = pd.DataFrame()
    for subdir, _, files in os.walk(input_dir):
        for file in files:
            if 'scaffold_info' in file:
                labels = file.replace('.IS_scaffold_info.tsv', '')
                df = pd.read_csv(os.path.join(subdir, file), sep="\t")
                df["experiment"] = labels
                scaf = pd.concat([scaf, df], ignore_index=True)
    return scaf

def load_genome_info(input_dir):
    genomes = pd.DataFrame()
    for subdir, _, files in os.walk(input_dir):
        for file in files:
            if 'genome_info' in file:
                labels = file.replace('.IS_genome_info.tsv', '')
                df = pd.read_csv(os.path.join(subdir, file), sep="\t")
                df["experiment"] = labels
                genomes = pd.concat([genomes, df], ignore_index=True)
    genomes["genome"] = genomes['genome'].str.removesuffix('.fa')
    return genomes

def process_snv_data(snv, scaf, stb_file, cutting_edge, coverage_limit, ratio_reads):
    stb = pd.read_table(stb_file, header=None, names=["scaffold", "genome"], index_col="scaffold")
    stb['genome'] = stb['genome'].str.removesuffix('.fa')
    snv['genome'] = stb.loc[snv['scaffold']].get('genome').values
    
    snv = snv[(snv.ref_base != 'N')]  # Remove unknown bases
    f = pd.merge(snv, scaf, on=["scaffold", "experiment"])
    
    f["ratio_VR"] = f.apply(lambda row: 1 if getattr(row, row.var_base) == 0 else getattr(row, row.var_base) / getattr(row, row.ref_base), axis=1)
    f["diff_cov"] = f["coverage"] - f["position_coverage"]
    f["end_limit"] = f["length"] - cutting_edge
    
    final_SNV = f[(f["position"] > cutting_edge) & (f["position"] < f["end_limit"]) & (f["diff_cov"].between(-coverage_limit, coverage_limit))]
    trustable = final_SNV[final_SNV["ratio_VR"] >= ratio_reads]["id"].unique()
    final_SNV = final_SNV[final_SNV['id'].isin(trustable)]
    
    return final_SNV

def main():
    parser = argparse.ArgumentParser(description="Process SNV data from multiple TSV files.")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing TSV files")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for filtered SNV file")
    parser.add_argument("-s", "--scaffold_to_bin", required=True, help="Scaffold to bin mapping file")
    parser.add_argument("-c", "--cutting_edges", type=int, default=100, help="Number of bases from the beginning and end of each scaffold to ignore")
    parser.add_argument("-l", "--coverage_limit", type=int, default=200, help="Maximum coverage difference between between scaffold and SNV position")
    parser.add_argument("-r", "--ratio_reads", type=float, default=0.10, help="Minimum ratio of variant base over reference base reads in order to considere a SNV high-confidence")
    args = parser.parse_args()
    
    snv = load_snv_files(args.input_dir)
    scaf = load_scaffold_info(args.input_dir)
    genomes = load_genome_info(args.input_dir)
    
    final_SNV = process_snv_data(snv, scaf, args.scaffold_to_bin, args.cutting_edge, args.coverage_limit, args.ratio_reads)
    output_file = os.path.join(args.output_dir, "filtered_SNVs.tsv")
    final_SNV.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered SNVs saved to {output_file}")
    
if __name__ == "__main__":
    main()
