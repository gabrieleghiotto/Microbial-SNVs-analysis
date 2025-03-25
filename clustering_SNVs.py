import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import numpy as np

def filter_snv_data(filtered_SNV, snv_type):
    if snv_type == "nonsyn":
        filtered_SNV = filtered_SNV[filtered_SNV["mutation_type"] == "N"]
    elif snv_type == "syn":
        filtered_SNV = filtered_SNV[filtered_SNV["mutation_type"] == "S"]
    return filtered_SNV

def process_snv(filtered_SNV, output_dir, position_coverage_threshold, min_freq, max_freq, apply_freq_filter, color_threshold):
    os.makedirs(output_dir, exist_ok=True)
    MAGs_name = filtered_SNV["genome"].unique()
    
    for genome in MAGs_name:
        genome_dir = os.path.join(output_dir, genome)
        os.makedirs(genome_dir, exist_ok=True)
        
        m = filtered_SNV[filtered_SNV["genome"] == genome]
        m = m[m["position_coverage"] >= position_coverage_threshold]
        
        evolution = m.pivot(index="id", columns="experiment", values="frequency").fillna(0)
        df = evolution.copy()
        df = df.astype(float).round(3)
        
        df = df.loc[(df != 1).any(axis=1)]  # Remove SNVs that have frequency of 1 in all time points
        df = df.loc[~(df == 0).all(axis=1)]  # Remove SNVs that have frequency of 0 in all time points
        
        # Apply frequency filtering if enabled
        if apply_freq_filter:
            df = df[(df >= min_freq).all(axis=1) & (df <= max_freq).all(axis=1)]
        
        output_table = os.path.join(genome_dir, f"{genome}_all.tsv")
        df.to_csv(output_table, sep="\t")
        
        if len(df.index) >= 0:
            fig, ax = plt.subplots(figsize=(7, 5))
            plot1 = os.path.join(genome_dir, f"clustering_{genome}.png")
            Z = sch.linkage(df, method="ward", optimal_ordering=False)
            dendrogram = sch.dendrogram(Z, labels=df.index, color_threshold=color_threshold, no_labels=True)
            plt.axhline(y=color_threshold, c='grey', lw=1, linestyle='dashed')
            plt.ylabel("Euclidean distance")
            plt.savefig(plot1, dpi=600, bbox_inches='tight')
            plt.close()
            
            # Assign unique strain labels based on color groups
            unique_colors = list(set(dendrogram["leaves_color_list"]))
            color_to_strain = {color: f"str{i}" for i, color in enumerate(unique_colors)}
            strains = [color_to_strain[x] for x in dendrogram["leaves_color_list"]]
            labels = dendrogram["ivl"]
            
            ceppo = pd.DataFrame(list(zip(strains, labels)), columns=['Strain', 'id']).set_index("id")
            melted_df = df.reset_index().melt(id_vars='id', var_name='timepoint', value_name='values').set_index('id')
            melted_df = pd.merge(melted_df, ceppo, left_index=True, right_index=True)
            
            fig, ax = plt.subplots(figsize=(14, 5))
            sns.lineplot(data=melted_df, x="timepoint", y="values", units="id", estimator=None, sort=False, legend=False, lw=.2, alpha=.2, hue="Strain", ax=ax)
            ax.set_ylim(0,1)
            plt.axhline(y=0.5, color='grey', linestyle='dashed')
            plt.xlabel("Time points")
            plt.ylabel("SNV frequency")
            plot2 = os.path.join(genome_dir, f"SNV_frequency_{genome}.png")
            plt.savefig(plot2, dpi=600, bbox_inches='tight')
            plt.close()

def main():
    parser = argparse.ArgumentParser(description="Process filtered SNV data and generate outputs.")
    parser.add_argument("-i", "--input_file", required=True, help="Input filtered SNV file")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for results")
    parser.add_argument("-t", "--snv_type", choices=["all", "syn", "nonsyn"], default="all", help="Type of SNVs to analyze")
    parser.add_argument("-p", "--position_coverage", type=int, default=10, help="Minimum position coverage threshold")
    parser.add_argument("--min_freq", type=float, default=0.05, help="Minimum SNV frequency threshold")
    parser.add_argument("--max_freq", type=float, default=0.95, help="Maximum SNV frequency threshold")
    parser.add_argument("--apply_freq_filter", action="store_true", help="Enable frequency filtering")
    parser.add_argument("--color_threshold", type=float, default=20, help="Color threshold for dendrogram clustering")
    args = parser.parse_args()
    
    filtered_SNV = pd.read_csv(args.input_file, sep="\t")
    filtered_SNV = filter_snv_data(filtered_SNV, args.snv_type)
    process_snv(filtered_SNV, args.output_dir, args.position_coverage, args.min_freq, args.max_freq, args.apply_freq_filter, args.color_threshold)
    print("Processing complete.")

if __name__ == "__main__":
    main()


