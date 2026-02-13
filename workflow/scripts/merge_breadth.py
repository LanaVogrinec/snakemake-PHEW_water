#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import os

def parse_args():
    p = argparse.ArgumentParser(description="Merge per-sample breadth files into a single matrix")
    p.add_argument("--breadth-files", nargs="+", required=True, help="List of per-sample breadth files")
    p.add_argument("--out", required=True, help="Output breadth matrix file")
    return p.parse_args()

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def main():
    args = parse_args()

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
