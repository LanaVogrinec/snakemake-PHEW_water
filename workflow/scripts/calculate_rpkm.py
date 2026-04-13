#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import re
import os
from datetime import datetime


def parse_args():
    p = argparse.ArgumentParser(description="Calculate RPKM and presence/absence")
    p.add_argument("--input", required=True)
    p.add_argument("--output_rpkm", required=True)
    p.add_argument("--output_pa", required=True)
    return p.parse_args()


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


CONTIG_COL = "rep_contig"


def get_sample_cols(df):
    return [
        c for c in df.columns
        if c not in {
            "rep_contig",
            "rep_contig_length",
            "singleton",
            "contig_length"
        }
    ]


def extract_contig_length(contig_name):
    m = re.search(r"_length_(\d+)_cov_", contig_name)
    if m:
        return int(m.group(1))
    raise ValueError(f"Cannot parse contig length from '{contig_name}'")


def main():
    args = parse_args()

    log(f"Reading input: {args.input}")
    df = pd.read_csv(args.input, sep="\t")

    log("Extracting contig lengths")
    df["contig_length"] = df[CONTIG_COL].apply(extract_contig_length)

    log("Selecting sample columns")
    read_cols = get_sample_cols(df)

    if len(read_cols) == 0:
        raise ValueError("No sample columns detected")

    log("Calculating total reads per sample")
    total_reads = df[read_cols].sum(axis=0)

    log("Calculating RPKM")
    rpkm_df = df.copy()

    for col in read_cols:
        rpkm_df[col] = (
            df[col] / df["contig_length"] * 1e3 / (total_reads[col] / 1e6)
        )

    log("Reshaping to long format")
    rpkm_long = rpkm_df.melt(
        id_vars=[CONTIG_COL],
        value_vars=read_cols,
        var_name="sample_id",
        value_name="RPKM"
    )

    log("Creating presence/absence table")
    pa_long = rpkm_long.copy()
    pa_long["presence"] = (pa_long["RPKM"] > 0).astype(int)
    pa_long = pa_long[[CONTIG_COL, "sample_id", "presence"]]

    log("Writing outputs")
    os.makedirs(os.path.dirname(args.output_rpkm), exist_ok=True)

    rpkm_long.to_csv(args.output_rpkm, sep="\t", index=False)
    pa_long.to_csv(args.output_pa, sep="\t", index=False)

    log("Done")


if __name__ == "__main__":
    main()