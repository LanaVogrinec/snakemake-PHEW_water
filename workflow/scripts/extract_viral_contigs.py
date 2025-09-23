#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Extract contigs classified as Viruses from merged FASTA based on taxonomy.")
parser.add_argument("--fasta", required=True, help="Input merged FASTA file with renamed contigs")
parser.add_argument("--taxonomy", required=True, help="Merged taxonomy TSV file from mmseqs2")
parser.add_argument("--output", required=True, help="Output FASTA file with only viral contigs")
args = parser.parse_args()

# Step 1: Load the taxonomy file (allow rows with and without lineage)
df = pd.read_csv(args.taxonomy, sep="\t", header=None, dtype=str)

# Step 2: Keep only rows with 9 fields (i.e., lineage is present)
df = df[df.apply(lambda row: len(row) == 9, axis=1)].copy()
df.columns = ["contig", "tax_id", "rank", "name", "retained", "assigned", "label_match", "support", "lineage"]

# Step 3: Extract domain from lineage (cleaning -_ prefix)
df["domain"] = df["lineage"].str.split(";").str[0].str.replace(r"^[-_]+", "", regex=True)

# Step 4: Select contigs where domain is "Viruses"
viral_ids = set(df[df["domain"] == "Viruses"]["contig"])

print(f"Total contigs with lineage: {len(df)}")
print(f"Total contigs classified as 'Viruses': {len(viral_ids)}")

# Step 5: Parse merged FASTA and write only matching sequences
count_written = 0
with open(args.output, "w") as out_fasta:
    for record in SeqIO.parse(args.fasta, "fasta"):
        if record.id in viral_ids:
            SeqIO.write(record, out_fasta, "fasta")
            count_written += 1

print(f"Total contigs written to output: {count_written}")
