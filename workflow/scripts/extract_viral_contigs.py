#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
from datetime import datetime

def parse_args():
    p = argparse.ArgumentParser(description="Extract contigs classified as Viruses from merged FASTA based on taxonomy")
    p.add_argument("--fasta", required=True, help="Input merged FASTA file with renamed contigs")
    p.add_argument("--taxonomy", required=True, help="Merged taxonomy TSV file from mmseqs2")
    p.add_argument("--output", required=True, help="Output FASTA file with only viral contigs")
    return p.parse_args()

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def main():
    args = parse_args()
    log(f"Reading taxonomy table: {args.taxonomy}")
    df = pd.read_csv(args.taxonomy, sep="\t", header=None, dtype=str)

    df = df[df.apply(lambda row: len(row) == 9, axis=1)].copy()
    df.columns = ["contig", "tax_id", "rank", "name", "retained", "assigned", "label_match", "support", "lineage"]

    df["domain"] = df["lineage"].str.split(";").str[0].str.replace(r"^[-_]+", "", regex=True)
    viral_ids = set(df[df["domain"] == "Viruses"]["contig"])

    log(f"Contigs with lineage: {len(df)}, classified as 'Viruses': {len(viral_ids)}")

    count_written = 0
    with open(args.output, "w") as out_fasta:
        for record in SeqIO.parse(args.fasta, "fasta"):
            if record.id in viral_ids:
                SeqIO.write(record, out_fasta, "fasta")
                count_written += 1

    log(f"Total contigs written to output: {count_written}")

if __name__ == "__main__":
    main()
