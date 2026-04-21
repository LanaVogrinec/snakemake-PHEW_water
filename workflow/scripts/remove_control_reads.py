#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import argparse


def parse_args():
    p = argparse.ArgumentParser(description="Remove carry-over reads using rep_contig matching")
    p.add_argument("--samples", required=True)
    p.add_argument("--carryovers", required=True)
    p.add_argument("--mapping", required=True)
    p.add_argument("--out-filtered", required=True)
    return p.parse_args()


def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)


def ensure_int(df):
    return df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def apply_carryover(samples_df, carry_df, mapping_df):
    out = samples_df.copy()

    for _, row in mapping_df.iterrows():
        sample_col = row["sample"]
        carry_col = row["carry_over"]

        if sample_col in out.columns and carry_col in carry_df.columns:
            sub = carry_df[carry_col].reindex(out.index).fillna(0)
            out[sample_col] = (out[sample_col] - sub).clip(lower=0)

    return out


def main():
    args = parse_args()

    log("Reading samples table")
    samples = pd.read_csv(args.samples, sep="\t").set_index("rep_contig")

    log("Reading carry-over table")
    carry = pd.read_csv(args.carryovers, sep="\t").set_index("rep_contig")

    log("Reading mapping table")
    mapping = pd.read_csv(args.mapping, sep="\t")

    samples = ensure_int(samples)
    carry = ensure_int(carry)

    log("Applying carry-over subtraction")
    filtered = apply_carryover(samples, carry, mapping)

    filtered = filtered.reset_index()
    filtered.to_csv(args.out_filtered, sep="\t", index=False)

    log("Done.")


if __name__ == "__main__":
    main()