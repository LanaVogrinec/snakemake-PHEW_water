#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
import glob
import os
import argparse
import time

parser = argparse.ArgumentParser(description="Merge taxonomy with representatives and combine read counts.")
parser.add_argument("--viral-reps", required=True, help="Representative contigs TSV (minlen filtered).")
parser.add_argument("--taxonomy", required=True, help="Merged taxonomy file for all contigs.")
parser.add_argument("--mapping-root", required=True, help="Root directory containing per-sample mapping TSVs.")
parser.add_argument("--out-reps-tax", required=True, help="Output TSV with representatives + taxonomy.")
parser.add_argument("--out-samples", required=True, help="Output read counts for environmental samples.")
parser.add_argument("--out-NKIs", required=True, help="Output read counts for NKIs / negative controls.")
parser.add_argument("--out-carryovers", required=True, help="Output read counts for carry-over samples.")
args = parser.parse_args()

# -----------------------------
# Step 1: Merge taxonomy
# -----------------------------
print(f"Reading representative contigs: {args.viral_reps}")
reps_df = pd.read_csv(args.viral_reps, sep="\t")
print(f"Loaded {len(reps_df)} representative contigs")

print(f"Reading taxonomy file: {args.taxonomy}")
tax_df = pd.read_csv(args.taxonomy, sep="\t", header=None)
tax_df.rename(columns={0: "rep_contig"}, inplace=True)  # first column is contig ID

# Merge by rep_contig
reps_tax_df = reps_df.merge(tax_df, how="left", on="rep_contig")
print(f"Merged taxonomy: {len(reps_tax_df)} rows")

# Write merged reps + taxonomy
reps_tax_df.to_csv(args.out_reps_tax, sep="\t", index=False)
print(f"Wrote representatives with taxonomy: {args.out_reps_tax}")

# -----------------------------
# Step 2: Combine read counts per sample
# -----------------------------
# Prepare separate dataframes
viral_samples = reps_tax_df.copy()
viral_NKIs = reps_tax_df.copy()
viral_carryovers = reps_tax_df.copy()

# Find mapping files recursively
pattern = os.path.join(args.mapping_root, "**", "08_*_coverm_bbmap_rep_all_minlen_500_noApr24_ANI90_AF85.tsv")
mapping_files = sorted(glob.glob(pattern, recursive=True))
print(f"Found {len(mapping_files)} mapping files to process")

if not mapping_files:
    raise FileNotFoundError(f"No mapping files found with pattern: {pattern}")

# Process each mapping file
for i, mapping_file in enumerate(mapping_files, 1):
    t0 = time.time()
    fname = os.path.basename(mapping_file)

    # Extract sample name
    sample_match = re.search(r'08_(.+?)_coverm', fname)
    if not sample_match:
        print(f"[WARN] Could not extract sample name from {mapping_file}, skipping.")
        continue
    sample_name = sample_match.group(1)
    new_col = f"{sample_name}_read_count"
    print(f"[{i}/{len(mapping_files)}] Processing sample: {sample_name}")

    try:
        mapping = pd.read_csv(mapping_file, sep="\t", usecols=[0, 1])
        mapping.columns = ["rep_contig", new_col]

        # Decide category
        if re.match(r'P2[34]_W\d{3}_E', sample_name):
            df_target = viral_samples
        elif re.match(r'P2[34]_N?KI', sample_name):
            df_target = viral_NKIs
        elif re.match(r'P2[34]_W\d{3}_C', sample_name):
            df_target = viral_carryovers
        else:
            print(f"[WARN] Sample {sample_name} did not match any category, skipping.")
            continue

        df_target = df_target.merge(mapping, how="left", on="rep_contig")

        # Update the correct DataFrame
        if re.match(r'P2[34]_W\d{3}_E', sample_name):
            viral_samples = df_target
        elif re.match(r'P2[34]_N?KI', sample_name):
            viral_NKIs = df_target
        elif re.match(r'P2[34]_W\d{3}_C', sample_name):
            viral_carryovers = df_target

        print(f"    Added read counts for {sample_name} ({time.time() - t0:.1f}s)")

    except Exception as e:
        print(f"[ERROR] Failed processing {mapping_file}: {e}")

    print("-"*80)

# -----------------------------
# Write outputs
# -----------------------------
print(f"Writing output for samples: {args.out_samples}")
viral_samples.to_csv(args.out_samples, sep="\t", index=False)

print(f"Writing output for NKIs: {args.out_NKIs}")
viral_NKIs.to_csv(args.out_NKIs, sep="\t", index=False)

print(f"Writing output for carry-overs: {args.out_carryovers}")
viral_carryovers.to_csv(args.out_carryovers, sep="\t", index=False)

print("All read count tables created successfully.")
