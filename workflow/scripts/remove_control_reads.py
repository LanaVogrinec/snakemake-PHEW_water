#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import argparse

def parse_args():
    p = argparse.ArgumentParser(description="Filter vOTUs by carry-over table")
    p.add_argument("--samples", required=True, help="Samples read counts table")
    p.add_argument("--carryovers", required=True, help="Carry-over read counts table")
    p.add_argument("--mapping", required=True, help="Carry-over mapping table")
    p.add_argument("--out-filtered", required=True, help="Filtered output table")
    return p.parse_args()

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def get_readcount_cols(df):
    return [c for c in df.columns if c.endswith("_read_count")]

def ensure_int(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    return df

def apply_carryover(counts_df, carry_df, mapping_df):
    out = counts_df.copy()
    for _, row in mapping_df.iterrows():
        sample_col = row["sample"] + "_read_count"
        carry_col  = row["carry_over"] + "_read_count"
        if sample_col in out.columns and carry_col in carry_df.columns:
            sub = carry_df[carry_col].reindex(out.index).fillna(0).astype(int)
            out[sample_col] = (out[sample_col] - sub).clip(lower=0)
    return out

def main():
    args = parse_args()
    log("Reading samples table")
    samples = pd.read_csv(args.samples, sep="\t").set_index("cluster_id")
    sample_cols = get_readcount_cols(samples)
    samples = ensure_int(samples, sample_cols)

    log("Reading carry-over table and mapping")
    carry = pd.read_csv(args.carryovers, sep="\t").set_index("cluster_id")
    mapping = pd.read_csv(args.mapping, sep="\t")
    carry_cols = get_readcount_cols(carry)
    carry = ensure_int(carry, carry_cols)

    log("Applying carry-over filtering")
    filtered_counts = apply_carryover(samples[sample_cols], carry, mapping)

    filtered_out = samples.copy()
    filtered_out[sample_cols] = filtered_counts
    filtered_out.insert(0, "cluster_id", filtered_out.index)
    filtered_out.to_csv(args.out_filtered, sep="\t", index=False)
    log(f"Wrote filtered table: {args.out_filtered}")
    log("Done.")

if __name__ == "__main__":
    main()