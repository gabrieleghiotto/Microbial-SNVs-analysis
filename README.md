# Microbial-SNVs-analysis

The UHGV is a comprehensive genomic resource of viruses from the human microbiome. Genomes were derived from [12 independent data sources](#data-sources) and annotated using a [uniform bioinformatics pipeline](#bioinformatics-pipeline):

## Table of contents
1. [Methods](#methods)
   * [Dependencies](#Dependencies)
   * [Bioinformatics pipeline](#bioinformatics-pipeline)
2. [Data availability](#data-availability)
   * [Recommended files](#recommended-files)
   * [All available files](#all-available-files)
3. [Bioinformatics tools that use the UHGV](#code-availability) 
   * [Contig-level taxonomic classification](CLASSIFY.md)
   * [Read-level abundance profiling](#read-level-abundance-profiling-with-phanta)
   * [Genome visualization](#genome-visualization)
      

## Methods

### Dependencies

The scripts require, beside the output from InStrain in a unique folder, the following dependencies: 

- Python:(https://www.python.org/)
- Pandas:(https://pandas.pydata.org/)
- Seaborn:(https://seaborn.pydata.org/)
- Matplotlib:(https://matplotlib.org/)
- Scipy:(https://scipy.org/)
- Numpy:(https://numpy.org/)

### Bioinformatics pipeline to generate input data

Sequences from these studies were combined and run through the following bioinformatics pipeline:
- [Spades](https://github.com/ablab/spades) was used for metagenomic assembly
- [Metabat](https://bitbucket.org/berkeleylab/metabat/src/master/), [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/), [Maxbin2](https://sourceforge.net/projects/maxbin2/) and [Vamb](https://github.com/RasmussenLab/vamb) were used for binning
- [checkM2](https://github.com/chklovski/CheckM2) was used to estimate completeness and contamination of recovered bins
- [dRep](https://github.com/MrOlm/drep) was used to dereplicate bins into MAGs, filtering for only medium-to-high quality genomes according to MIMAG standards
- [GTDB r214](https://gtdb.ecogenomic.org/) and [GTDB-tk](https://github.com/Ecogenomics/GTDBTk) were used to assign taxonomy to all prokaryotic genomes
- [Prodigal](https://github.com/hyattpd/Prodigal) was used to identify protein-coding genes
- [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) was used for gene functional annotation
- [parse_stb.py](https://github.com/MrOlm/drep/blob/master/helper_scripts/) was used to generate a scaffold to bin file
- [InStrain](https://github.com/MrOlm/inStrain) was used to call SNVs on MAGs
- [Bowtie2](https://github.com/BenLangmead/bowtie2) was used to align short reads to the assembly and to the MAGs

For additional details, please refer to our manuscript: (in preparation).

## Data availability

The entire resource has been deposited at SRA () and is  available upon request.

### Recommended files
For most analyses of InStrain output, we recommend using these files:
- [scaffold_info.tsv]
- [SNVs.tsv]
- [genome_info.tsv]

### All available files:

- [linkage.tsv]
- [gene_info.tsv]

## Code availability
The repository is composed of two Python scripts at the moment:
- [filter_SNVs.py] takes as input files the output from InStrain profile, on one or multiple metagenomic samples
- [cluster_SNVs.py] takes the filtered dataset generated with [filter_SNV.py] and performs further filtering of the SNVs (if wanted) and then clustering of SNVs based on their frequency over time
