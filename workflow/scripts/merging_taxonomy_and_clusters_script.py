#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s:%(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

parser = argparse.ArgumentParser(description="Merge cluster assignments with merged MMseqs2 taxonomy data.")
parser.add_argument("--clusters", required=True, help="TSV file with contig_id and cluster_id.")
parser.add_argument("--taxonomy", required=True, help="Merged MMseqs2 taxonomy file.")
parser.add_argument("--output", required=True, help="Output merged TSV file.")
args = parser.parse_args()

# Load clusters
logging.info("Loading cluster file...")
clusters_df = pd.read_csv(args.clusters, sep="\t")
if not {"contig_id", "cluster_id"}.issubset(clusters_df.columns):
    raise ValueError("clusters file must contain columns: contig_id, cluster_id")
logging.info(f"Loaded {len(clusters_df)} cluster assignments.")

# Compute singleton flag
cluster_sizes = clusters_df.groupby("cluster_id")["contig_id"].size().rename("cluster_size")
clusters_df = clusters_df.merge(cluster_sizes, on="cluster_id", how="left")
clusters_df["singleton"] = (clusters_df["cluster_size"] == 1).astype(int)

# Load merged taxonomy
logging.info(f"Reading taxonomy file: {args.taxonomy}")
taxonomy_df = pd.read_csv(args.taxonomy, sep="\t", header=None)
taxonomy_df.columns = [
    "contig_id", "taxonomic_identifier", "taxonomic_rank", "taxonomic_name",
    "no_of_retained_fragments", "no_of_assigned_fragments",
    "in_agreement_with_contig_label", "support", "lineage"
]
logging.info(f"Total taxonomy entries loaded: {len(taxonomy_df)}")

# Merge
merged_df = clusters_df.merge(taxonomy_df, on="contig_id", how="left")
logging.info(f"Merged table has {len(merged_df)} rows.")

# Optional: drop cluster_size if you only want singleton 1/0
merged_df = merged_df.drop(columns=["cluster_size"])

# Write
merged_df.to_csv(args.output, sep="\t", index=False)
logging.info(f"Merged table written to {args.output}")
logging.info("Merge complete.")
