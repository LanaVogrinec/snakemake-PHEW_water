#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args():
    p = argparse.ArgumentParser(description="Select longest representative contig per cluster")
    p.add_argument("--fasta", required=True)
    p.add_argument("--clusters", required=True)
    p.add_argument("--reps-fasta", required=True)
    p.add_argument("--reps-tsv", required=True)
    p.add_argument("--stats", required=True)
    return p.parse_args()


def log(msg):
    print(f"[INFO] {msg}", flush=True)


def main():
    args = parse_args()

    log("Reading cluster table")
    df = pd.read_csv(args.clusters, sep="\t")

    if not {"contig_id", "cluster_id"}.issubset(df.columns):
        raise ValueError("Missing contig_id or cluster_id")

    log("Loading FASTA")
    seqs = {r.id: r for r in SeqIO.parse(args.fasta, "fasta")}

    df["length"] = df["contig_id"].map(lambda x: len(seqs[x].seq) if x in seqs else 0)
    df = df[df["length"] > 0].copy()

    log("Selecting longest contig per cluster")

    # keep all columns explicitly so nothing is dropped
    df = df.sort_values("length", ascending=False)
    reps = df.loc[df.groupby("cluster_id")["length"].idxmax()].copy()

    # build final representative ID
    reps["rep_contig"] = reps["contig_id"] + "_" + reps["cluster_id"]

    # reorder columns
    reps_df = reps[["rep_contig", "contig_id", "cluster_id", "length"]].rename(
        columns={"length": "rep_contig_length"}
    )

    log("Writing TSV")
    reps_df.to_csv(args.reps_tsv, sep="\t", index=False)

    log("Writing FASTA")
    with open(args.reps_fasta, "w") as out:
        for _, row in reps_df.iterrows():
            rec = seqs[row["contig_id"]]
            new_id = row["rep_contig"]
            SeqIO.write(
                SeqRecord(rec.seq, id=new_id, name=new_id, description=""),
                out,
                "fasta"
            )

    log("Writing stats")
    with open(args.stats, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_clusters\t{reps_df['cluster_id'].nunique()}\n")


if __name__ == "__main__":
    main()