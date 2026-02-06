#!/usr/bin/env python3

import argparse
import logging
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def main():
    parser = argparse.ArgumentParser(
        description="Select longest representative contig per cluster and write FASTA + TSV + stats."
    )
    parser.add_argument("--fasta", required=True,
                        help="FASTA file with all viral contigs.")
    parser.add_argument("--clusters", required=True,
                        help="TSV file with columns: contig_id, cluster_id.")
    parser.add_argument("--reps-fasta", required=True,
                        help="Output FASTA file with representative contigs.")
    parser.add_argument("--reps-tsv", required=True,
                        help="Output TSV file with representative metadata.")
    parser.add_argument("--stats", required=True,
                        help="Output TSV/TXT file with cluster statistics.")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load clusters table
    # ------------------------------------------------------------------
    clusters_df = pd.read_csv(args.clusters, sep="\t")
    if not {"contig_id", "cluster_id"}.issubset(clusters_df.columns):
        raise ValueError("Clusters file must contain columns: contig_id, cluster_id")

    logging.info(f"Loaded {len(clusters_df)} contig cluster assignments")

    # Cluster sizes
    cluster_sizes = (
        clusters_df
        .groupby("cluster_id")
        .size()
        .rename("cluster_size")
        .reset_index()
    )

    # ------------------------------------------------------------------
    # Load FASTA
    # ------------------------------------------------------------------
    seq_dict = {rec.id: rec for rec in SeqIO.parse(args.fasta, "fasta")}
    logging.info(f"Loaded {len(seq_dict)} contigs from FASTA")

    # Warn about missing contigs
    missing = set(clusters_df["contig_id"]) - set(seq_dict.keys())
    if missing:
        logging.warning(f"{len(missing)} contigs from clusters file not found in FASTA")

    # ------------------------------------------------------------------
    # Add contig length
    # ------------------------------------------------------------------
    clusters_df["length"] = clusters_df["contig_id"].map(
        lambda cid: len(seq_dict[cid].seq) if cid in seq_dict else 0
    )

    before = len(clusters_df)
    clusters_df = clusters_df[clusters_df["length"] > 0].copy()
    after = len(clusters_df)
    if after < before:
        logging.warning(f"Removed {before - after} contigs with length 0 (missing in FASTA)")

    # ------------------------------------------------------------------
    # Pick longest contig per cluster
    # ------------------------------------------------------------------
    reps_df = (
        clusters_df
        .sort_values("length", ascending=False)
        .groupby("cluster_id", as_index=False)
        .first()
    )

    reps_df = reps_df[["cluster_id", "contig_id", "length"]]
    reps_df = reps_df.rename(columns={
        "contig_id": "rep_contig",
        "length": "rep_contig_length"
    })

    # ------------------------------------------------------------------
    # Add singleton annotation
    # ------------------------------------------------------------------
    reps_df = reps_df.merge(cluster_sizes, on="cluster_id", how="left")
    reps_df["singleton"] = (reps_df["cluster_size"] == 1).astype(int)
    reps_df = reps_df.drop(columns=["cluster_size"])

    logging.info(f"Selected {len(reps_df)} representative contigs (one per cluster)")

    # ------------------------------------------------------------------
    # Write representatives TSV
    # ------------------------------------------------------------------
    reps_df.to_csv(args.reps_tsv, sep="\t", index=False)
    logging.info(f"Wrote representatives table to {args.reps_tsv}")

    # ------------------------------------------------------------------
    # Write FASTA
    # ------------------------------------------------------------------
    with open(args.reps_fasta, "w") as out_fa:
        for _, row in reps_df.iterrows():
            cluster_id = row["cluster_id"]
            contig_id = row["rep_contig"]
            rec = seq_dict[contig_id]
            new_id = f"{contig_id}_{cluster_id}"

            new_rec = SeqRecord(
                rec.seq,
                id=new_id,
                name=new_id,
                description=""
            )
            SeqIO.write(new_rec, out_fa, "fasta")

    logging.info(f"Wrote representative contigs FASTA to {args.reps_fasta}")

    # ------------------------------------------------------------------
    # Write stats file
    # ------------------------------------------------------------------
    total_clusters = len(reps_df)
    n_singletons = int((reps_df["singleton"] == 1).sum())
    n_non_singletons = total_clusters - n_singletons

    with open(args.stats, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"total_clusters\t{total_clusters}\n")
        fh.write(f"singletons\t{n_singletons}\n")
        fh.write(f"non_singletons\t{n_non_singletons}\n")

    logging.info(f"Wrote cluster statistics to {args.stats}")


if __name__ == "__main__":
    main()
