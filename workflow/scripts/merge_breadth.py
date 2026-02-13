#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import os

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def main():
    parser = argparse.ArgumentParser(description="Merge per-sample breadth files into a single matrix")
    parser.add_argument("--breadth-files", nargs="+", required=True, help="List of per-sample breadth files")
    parser.add_argument("--out", required=True, help="Output breadth matrix file")
    args = parser.parse_args()

    if len(args.breadth_files) == 0:
        raise FileNotFoundError("No breadth files provided!")

    log(f"Found {len(args.breadth_files)} breadth files")

    dfs = []
    for f in args.breadth_files:
        sample = os.path.basename(os.path.dirname(f))
        df = pd.read_csv(f, sep="\t", header=None, names=["contig", sample])
        dfs.append(df)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="contig", how="outer")

    merged.fillna(0, inplace=True)
    merged.to_csv(args.out, sep="\t", index=False)
    log(f"Finished: {args.out}")

if __name__ == "__main__":
    main()
