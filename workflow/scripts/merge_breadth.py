#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Merge per-sample breadth files into a single matrix.")
    parser.add_argument(
        "--breadth-files",
        nargs="+",
        required=True,
        help="List of per-sample breadth TSV files."
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file for merged breadth matrix."
    )
    return parser.parse_args()

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def main():
    args = parse_args()

    dfs = []
    for f in args.breadth_files:
        # Extract sample name from two levels up (folder containing test_Apr24_DS1)
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(f)))
        log(f"Reading {f} as sample '{sample_name}'")
        df = pd.read_csv(f, sep="\t", header=None, names=["contig", sample_name])
        df = df.set_index("contig")
        dfs.append(df)

    log(f"Merging {len(dfs)} breadth files")
    merged = pd.concat(dfs, axis=1)
    merged = merged.fillna(0)

    log(f"Writing merged breadth matrix to {args.out}")
    merged.to_csv(args.out, sep="\t")

    log("Done!")

if __name__ == "__main__":
    main()