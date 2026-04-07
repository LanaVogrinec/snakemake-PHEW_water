#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import argparse

def parse_args():
    p = argparse.ArgumentParser(description="Filter vOTUs by min reads and breadth")
    p.add_argument("--samples", required=True, help="Samples read counts table (carry-over cleaned)")
    p.add_argument("--breadth", required=True, help="Breadth table")
    p.add_argument("--out-filtered", required=True, help="Filtered output table")
    p.add_argument("--out-stats", required=True, help="Stats output table")
    p.add_argument("--min-reads", type=int, default=2)
    p.add_argument("--min-breadth", type=float, default=0.0)
    return p.parse_args()

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def get_readcount_cols(df):
    return [c for c in df.columns if c.endswith("_read_count")]

def ensure_int(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    return df

def apply_read_breadth_filter(counts_df, rep_contig_series, breadth_df, min_reads, min_breadth):
    out = counts_df.copy()
    for sample in out.columns:
        sample_id = sample.replace("_read_count", "")
        if sample_id not in breadth_df.columns:
            continue
        breadth_vals = rep_contig_series.map(breadth_df[sample_id]).fillna(0)
        low_reads = out[sample] < min_reads
        low_breadth = breadth_vals < min_breadth
        out.loc[low_reads | low_breadth, sample] = 0
    return out

def count_vOTU_detections(df):
    return (df > 0).sum().sum()

def main():
    args = parse_args()
    log("Reading samples table")
    samples = pd.read_csv(args.samples, sep="\t").set_index("cluster_id")
    sample_cols = get_readcount_cols(samples)
    samples = ensure_int(samples, sample_cols)

    log("Reading breadth table")
    breadth = pd.read_csv(args.breadth, sep="\t")
    breadth = breadth[breadth['contig'] != 'contig'].copy()
    breadth['contig'] = breadth['contig'].astype(str).str.replace(r'_cluster_\d+$', '', regex=True)
    breadth = breadth.set_index('contig')
    breadth = breadth.apply(pd.to_numeric, errors='coerce').fillna(0)

    rep_contigs = samples["rep_contig"].astype(str).str.replace(r'_cluster_\d+$', '', regex=True)

    log(f"Applying min_reads={args.min_reads} and min_breadth={args.min_breadth} filtering")
    filtered_counts = apply_read_breadth_filter(samples[sample_cols], rep_contigs, breadth,
                                                args.min_reads, args.min_breadth)

    filtered_out = samples.copy()
    filtered_out[sample_cols] = filtered_counts
    filtered_out.insert(0, "cluster_id", filtered_out.index)
    filtered_out.to_csv(args.out_filtered, sep="\t", index=False)
    log(f"Wrote filtered table: {args.out_filtered}")

    stats = pd.DataFrame({
        "step": ["after_read_breadth"],
        "vOTU_detections": [count_vOTU_detections(filtered_counts)]
    })
    stats.to_csv(args.out_stats, sep="\t", index=False)
    log(f"Wrote stats table: {args.out_stats}")
    log("Done.")

if __name__ == "__main__":
    main()