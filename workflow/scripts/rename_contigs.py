#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO

def parse_args():
    p = argparse.ArgumentParser(description="Rename contigs in a FASTA by prepending a prefix")
    p.add_argument("--input", required=True, help="Input FASTA file")
    p.add_argument("--output", required=True, help="Output FASTA file with renamed headers")
    p.add_argument("--prefix", required=True, help="Prefix to prepend to each contig name")
    return p.parse_args()

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def main():
    args = parse_args()
    log(f"Reading {args.input} and writing renamed contigs to {args.output}")

    with open(args.input) as infile, open(args.output, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            record.id = f"{args.prefix}_{record.id}"
            record.description = ""
            SeqIO.write(record, outfile, "fasta")

    log("Done.")

if __name__ == "__main__":
    main()
