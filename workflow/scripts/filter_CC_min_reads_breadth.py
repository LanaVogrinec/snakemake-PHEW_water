#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import argparse

# ============================================================
# HELPERS
# ============================================================
def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def get_readcount_cols(df):
    return [c for c in df.columns if c.endswith("_read_count")]

def ensure_int(df, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    return df

def apply_carryover(counts_df, carry_df, mapping_df=None):
    """Subtract carryover reads per sample."""
    out = counts_df.copy()
    # mapping_df is ignored because carryovers already have same columns as counts
    for c in counts_df.columns:
        if c in carry_df.columns:
            out[c] = (out[c] - carry_df[c]).clip(lower=0)
    return out

def apply_read_breadth_filter(counts_df, rep_contig_series, breadth_df, min_reads, min_breadth):
    out = counts_df.copy()
    for sample in out.columns:
        if sample.replace("_read_count", "") not in breadth_df.columns:
            continue
        breadth_vals = rep_contig_series.map(breadth_df[sample.replace("_read_count", "")]).fillna(0)
        low_reads = out[sample] < min_reads
        low_breadth = breadth_vals < min_breadth
        out.loc[low_reads | low_breadth, sample] = 0
    return out

def count_vOTUs(df):
    """Count number of vOTUs with >=1 read summed across all samples."""
    return (df > 0).any(axis=1).sum()

# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="Filter vOTUs with carry-over and read/breadth thresholds")
    parser.add_argument("--samples", required=True)
    parser.add_argument("--carryovers", required=True)
    parser.add_argument("--breadth", required=True)
    parser.add_argument("--out-filtered", required=True)
    parser.add_argument("--out-stats", required=True)
    parser.add_argument("--min-reads", type=int, default=2)
    parser.add_argument("--min-breadth", type=float, default=0.0)
    args = parser.parse_args()

    log("Reading input tables")
    samples = pd.read_csv(args.samples, sep="\t")
    carry   = pd.read_csv(args.carryovers, sep="\t")
    breadth = pd.read_csv(args.breadth, sep="\t").set_index(breadth.columns[0])

    rep_contigs = samples.set_index("cluster_id")["rep_contig"]
    sample_cols = get_readcount_cols(samples)
    samples_indexed = samples.set_index("cluster_id")
    samples_indexed = ensure_int(samples_indexed, sample_cols)
    carry_indexed   = carry.set_index("cluster_id")
    carry_indexed   = ensure_int(carry_indexed, sample_cols)

    base_counts = samples_indexed[sample_cols]

    # ============================================================
    # Step 1: Carry-over filtering
    # ============================================================
    log("Applying carry-over filtering")
    cc_counts = apply_carryover(base_counts, carry_indexed)
    
    # ============================================================
    # Step 2: Min read/breadth filtering
    # ============================================================
    log(f"Applying min_reads={args.min_reads} and min_breadth={args.min_breadth} filtering")
    filtered_counts = apply_read_breadth_filter(cc_counts, rep_contigs, breadth, args.min_reads, args.min_breadth)

    # ============================================================
    # Save filtered table
    # ============================================================
    filtered_out = samples_indexed.copy()
    filtered_out[sample_cols] = filtered_counts[sample_cols]
    filtered_out.insert(0, "cluster_id", filtered_out.index)
    filtered_out.to_csv(args.out_filtered, sep="\t", index=False)
    log(f"Wrote filtered table: {args.out_filtered}")

    # ============================================================
    # Save stats
    # ============================================================
    stats = pd.DataFrame({
        "step": ["before_filtering", "after_carryover", "after_read_breadth"],
        "vOTUs_detected": [
            count_vOTUs(base_counts),
            count_vOTUs(cc_counts),
            count_vOTUs(filtered_counts)
        ]
    })
    stats.to_csv(args.out_stats, sep="\t", index=False)
    log(f"Wrote stats table: {args.out_stats}")
    log("Done.")

if __name__ == "__main__":
    main()
