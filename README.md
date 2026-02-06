# Snakemake workflow: `project PHEW`

A Snakemake workflow for virus detection in metagenomic data from various samples (plants, river water, wastewater). This workflow is based on the workflow from [project-tobamo] (https://github.com/nezapajek/project-tobamo) created by Neža Pajek (https://github.com/nezapajek). 

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

## Usage

To run this Snakemake workflow, follow these steps:

1. **Prerequisites**: Ensure that you have the following software and dependencies installed:
   - Snakemake (version ≥6.3.0)
   
Snakemake is best installed via the Mamba package manager (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can be installed via Mambaforge. Given that Mamba is installed, run:

        mamba create -c conda-forge -c bioconda --name snakemake snakemake

For all following commands ensure that this environment is activated via:

        conda activate snakemake

2. **Clone Repository**: Clone this repository to your local machine.

First, create an appropriate project working directory on your system and enter it:

        mkdir -p path/to/project-workdir
        cd path/to/project-workdir

Then, clone the github repository:
   
        git clone https://github.com/LanaVogrinec/snakemake-marchantia

3. **Configure workflow** Modify the workflow to your needs.

Modify the configuration file (config/config.yaml) and the Snakefile (workflow/Snakefile) to specify input files, parameters, and any other settings specific to your project. 

4. **Run workflow** Run the workflow with the following command: 

        snakemake --use-conda --cores <num_cores> 

## Workflow
Per sample:
1. Trimming of illumina adapters + QC with trimmomatic (trimmomatic.smk)
2. Trimming of preamplification adapters with cutadapt (cutadapt.smk)
3. Aseembly with SPAdes (--metaspades.py) (spades.smk)
4. Determine taxonomic classification for each contig with mmseqs2 (mmseqs2.smk)
5. Transform the results of the mmseqs2 (database) into a tsv file (mmseqs2_tsv.smk)
6. Rename contigs to include the sample name at the beginning (rename_contigs.smk)
7. Merge all contigs from all samples together into one file (merge_contigs.smk)

Per merged contigs
8. Merge the name of each contig with its associated taxonomy into one combined file (merge_taxonomy_tables.tsv)
9. Extract only contigs that were classified as viral with mmseqs ('Viruses') (extract_viral_contigs.smk)
10. Perform an all-vs-all blastn on all viral contigs (blast.smk)
11. Cluster the contigs based on blastn all-vs-all results and a network graph (clustering.smk)
12. Merge the information about cluster membership with the taxonomy per contig (merge_taxonomy_clusters.smk)
13. Split the file of clusters with taxonomy into only clusters > 1 contig and singleton clusters (split_clusters.smk)
14. Merge the information on pairwise hits from blastn with cluster membership (annotate_clusters_blast.smk)
--- somewhere here we need to include automated hybrid cluster breaking ? ---
15. Map reads per sample to all merged contigs, 90 % identity threshold (viral + non-viral) (bbmap_map.smk)
16. Filter only for reads that mapped with 95 % identity, 95 % aligned fraction (coverm_quant_bbmap.smk)
17. Select representatives of clusters based on length (longest contig per cluster) (select_cluster_representatives.smk)
18. Extract only reads mapped to viral contigs from per-sample mapping bam files and output into F/R fastq read files (extract_viral_reads.smk)
19. Map only viral reads to representative contigs per sample (bbmap_rep.smk)
20. Filter only for reads that mapped with 95 % identity, 95 % aligned fraction (coverm_quant_bbmap_rep.smk)


If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of the original repository.
