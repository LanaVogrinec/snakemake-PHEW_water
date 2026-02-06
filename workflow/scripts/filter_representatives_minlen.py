#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import logging
import pandas as pd
from Bio import SeqIO

# -------------------------------------------------------------------
# Logging setup
# -------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

# -------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Filter representative contigs by minimum length using TSV metadata."
    )
    parser.add_argument("--reps-fasta", required=True, help="FASTA of representative contigs.")
    parser.add_argument("--reps-tsv", required=True, help="TSV list of representative contigs.")
    parser.add_argument("--out-fasta", required=True, help="Output FASTA of filtered representatives.")
    parser.add_argument("--out-tsv", required=True, help="Output TSV list of filtered representatives.")
    parser.add_argument("--stats", required=True, help="Output stats file.")
    parser.add_argument("--min-len", type=int, default=500, help="Minimum contig length to keep.")
    args = parser.parse_args()

    logging.info(f"Minimum length cutoff: {args.min_len} nt")

    # -------------------------------------------------------------------
    # Load TSV
    # -------------------------------------------------------------------
    reps_df = pd.read_csv(args.reps_tsv, sep="\t")

    id_col = "rep_contig"
    len_col = "rep_contig_length"

    missing_cols = {id_col, len_col} - set(reps_df.columns)
    if missing_cols:
        sys.exit(f"ERROR: Missing required column(s) in TSV: {', '.join(missing_cols)}")

    singleton_col = "singleton"
    has_singleton = singleton_col in reps_df.columns

    logging.info(f"Using contig ID column: {id_col}")
    logging.info(f"Using contig length column: {len_col}")
    if has_singleton:
        logging.info(f"Including singleton column: {singleton_col}")

    # -------------------------------------------------------------------
    # Filter TSV by length
    # -------------------------------------------------------------------
    keep_df = reps_df[reps_df[len_col] >= args.min_len].copy()
    keep_ids = set(keep_df[id_col].astype(str))

    logging.info(f"Representatives in TSV: {len(reps_df)}")
    logging.info(f"Representatives kept after length filter: {len(keep_df)}")

    if len(keep_df) == 0:
        logging.warning(
            "No representatives passed the length filter. "
            "Check input TSV and --min-len value."
        )

    # Write filtered TSV (retain singleton column if present)
    keep_df.to_csv(args.out_tsv, sep="\t", index=False)
    logging.info(f"Wrote filtered representative list: {args.out_tsv}")

    # -------------------------------------------------------------------
    # Filter FASTA to match TSV
    # -------------------------------------------------------------------
    total_records = 0
    kept_records = 0
    with open(args.out_fasta, "w") as out_fa:
        for rec in SeqIO.parse(args.reps_fasta, "fasta"):
            total_records += 1
            # Strip _cluster_XX suffix to match TSV IDs
            rec_id_base = rec.id.split("_cluster_")[0]
            if rec_id_base in keep_ids:
                SeqIO.write(rec, out_fa, "fasta")
                kept_records += 1

    logging.info(f"FASTA records read: {total_records}")
    logging.info(f"FASTA records written: {kept_records}")
    logging.info(f"Wrote filtered representative FASTA: {args.out_fasta}")

    if kept_records != len(keep_ids):
        logging.warning(
            f"Mismatch: TSV kept {len(keep_ids)} IDs but FASTA wrote {kept_records}. "
            "Some representative IDs may be missing from the FASTA."
        )

    # -------------------------------------------------------------------
    # Write stats file
    # -------------------------------------------------------------------
    total_reps = len(keep_df)
    n_singletons = int((keep_df[singleton_col] == 1).sum()) if has_singleton else 0
    n_non_singletons = total_reps - n_singletons

    with open(args.stats, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"total_representatives\t{total_reps}\n")
        fh.write(f"singletons\t{n_singletons}\n")
        fh.write(f"non_singletons\t{n_non_singletons}\n")

    logging.info(f"Wrote filtered representatives stats to {args.stats}")
    logging.info("Done.")


# -------------------------------------------------------------------
# Entry point
# -------------------------------------------------------------------
if __name__ == "__main__":
    main()
