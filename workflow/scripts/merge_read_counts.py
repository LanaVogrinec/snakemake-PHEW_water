#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import re
import glob
import os
import sys
import time

# ------------------------
# USER-ADJUSTABLE PATTERNS
# ------------------------
SAMPLE_PATTERN = r'P2[34]_W\d{3}_E(_\d+)?$'
NKI_PATTERN    = r'P2[34]_N?KI(_\d+)?$'
CARRY_PATTERN  = r'P2[34]_W\d{3}_C(_\d+)?$'
# ------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Merge taxonomy with representatives and combine read counts.")
    p.add_argument("--viral-reps", required=True, help="Representative contigs TSV (minlen filtered).")
    p.add_argument("--taxonomy", required=True, help="Merged taxonomy file for all contigs.")
    p.add_argument("--mapping-root", required=True, help="Root directory containing per-sample mapping TSVs.")
    p.add_argument("--out-reps-tax", default=None, help="Output TSV with representatives + taxonomy.")
    p.add_argument("--out-samples", default=None, help="Output read counts for environmental samples.")
    p.add_argument("--out-NKIs", default=None, help="Output read counts for NKIs / negative controls.")
    p.add_argument("--out-carryovers", default=None, help="Output read counts for carry-over samples.")
    return p.parse_args()

def log(msg):
    print(f"[INFO] {msg}", flush=True)

def read_tsv_file(path, **kwargs):
    """Read TSV robustly with checks and Python engine."""
    if not os.path.exists(path):
        sys.exit(f"[ERROR] File does not exist: {path}")
    if os.path.getsize(path) == 0:
        sys.exit(f"[ERROR] File is empty: {path}")
    try:
        df = pd.read_csv(path, **kwargs, engine="python")
        return df
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read {path}: {e}")

def categorize_sample(sample_name):
    """Return category based on sample naming pattern"""
    if re.match(SAMPLE_PATTERN, sample_name):
        return "sample"
    elif re.match(NKI_PATTERN, sample_name):
        return "NKI"
    elif re.match(CARRY_PATTERN, sample_name):
        return "carryover"
    else:
        return None

def main():
    args = parse_args()

    log(f"Reading representative contigs: {args.viral_reps}")
    reps_df = read_tsv_file(args.viral_reps, sep="\t")
    log(f"Loaded {len(reps_df)} representative contigs")

    log(f"Reading taxonomy file: {args.taxonomy}")
    tax_df = read_tsv_file(args.taxonomy, sep="\t", header=None)
    tax_df.rename(columns={0: "rep_contig"}, inplace=True)
    log(f"Loaded {len(tax_df)} taxonomy entries")

    reps_tax_df = reps_df.merge(tax_df, how="left", on="rep_contig")

    if args.out_reps_tax:
        reps_tax_df.to_csv(args.out_reps_tax, sep="\t", index=False)
        log(f"Wrote representatives with taxonomy: {args.out_reps_tax}")

    viral_samples = reps_tax_df.copy()
    viral_NKIs = reps_tax_df.copy()
    viral_carryovers = reps_tax_df.copy()

    # Discover mapping files
    pattern = os.path.join(args.mapping_root, "**", "08_*_coverm_filtered_reps.tsv")
    mapping_files = sorted(glob.glob(pattern, recursive=True))
    log(f"Found {len(mapping_files)} mapping files")

    if not mapping_files:
        log(f"[WARN] No mapping files found with pattern: {pattern}. Exiting.")
        return

    for i, mapping_file in enumerate(mapping_files, 1):
        t0 = time.time()
        fname = os.path.basename(mapping_file)
        sample_match = re.search(r'08_(.+?)_coverm', fname)
        if not sample_match:
            log(f"[WARN] Could not extract sample name from {mapping_file}, skipping.")
            continue
        sample_name = sample_match.group(1)
        category = categorize_sample(sample_name)
        if category is None:
            log(f"[WARN] Sample {sample_name} did not match any category, skipping.")
            continue

        new_col = f"{sample_name}_read_count"
        log(f"[{i}/{len(mapping_files)}] Processing {category}: {sample_name}")

        try:
            mapping = read_tsv_file(mapping_file, sep="\t", usecols=[0,1])
            mapping.columns = ["Contig", new_col]
            mapping["Contig"] = mapping["Contig"].str.replace(r"_cluster_.*$", "", regex=True)

            if category == "sample":
                df_target = viral_samples
            elif category == "NKI":
                df_target = viral_NKIs
            elif category == "carryover":
                df_target = viral_carryovers

            df_target = df_target.merge(mapping, how="left", left_on="rep_contig", right_on="Contig")
            df_target.drop(columns=["Contig"], inplace=True)

            if category == "sample":
                viral_samples = df_target
            elif category == "NKI":
                viral_NKIs = df_target
            elif category == "carryover":
                viral_carryovers = df_target

            log(f"    Added read counts for {sample_name} ({time.time() - t0:.1f}s)")

        except Exception as e:
            log(f"[ERROR] Failed processing {mapping_file}: {e}")

    if args.out_samples:
        viral_samples.to_csv(args.out_samples, sep="\t", index=False)
        log(f"Wrote environmental samples read counts: {args.out_samples}")

    if args.out_NKIs:
        viral_NKIs.to_csv(args.out_NKIs, sep="\t", index=False)
        log(f"Wrote NKIs read counts: {args.out_NKIs}")

    if args.out_carryovers:
        viral_carryovers.to_csv(args.out_carryovers, sep="\t", index=False)
        log(f"Wrote carry-over read counts: {args.out_carryovers}")

    log("All requested read count tables created successfully.")

if __name__ == "__main__":
    main()