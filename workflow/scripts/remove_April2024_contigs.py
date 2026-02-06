#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip

# -----------------------------
# INPUT / OUTPUT (edit these)
# -----------------------------
IN_BLAST_GZ  = "/DATB/scratch/lanav/PHEW/snakemake-PHEW_water/results/merged/blastn_pairs_sorted_no-self.tsv.gz"
OUT_BLAST_GZ = "/DATB/scratch/lanav/PHEW/snakemake-PHEW_water/results/merged/blastn_pairs_sorted_no-self_no_april2024.tsv.gz"

# Samples to exclude (April 2024)
EXCLUDE_SAMPLES = {"P24_W001_E", "P24_W002_E", "P24_W003_E", "P24_W004_E"}


def sample_from_contig(contig_id: str) -> str:
    """
    Extract sample ID from contig ID.
    Assumes contigs begin with SAMPLEID, e.g. P24_W001_E_NODE_...
    """
    if "_NODE_" in contig_id:
        return contig_id.split("_NODE_", 1)[0]
    # fallback: take first 3 underscore-delimited parts (P24_W001_E)
    parts = contig_id.split("_")
    if len(parts) >= 3:
        return "_".join(parts[:3])
    return contig_id


def main():
    kept = 0
    skipped = 0

    with gzip.open(IN_BLAST_GZ, "rt") as fin, gzip.open(OUT_BLAST_GZ, "wt") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            cols = line.split("\t")

            # Expected format from your sorter:
            # A  B  qseqid  sseqid  pident  length  qstart  qend  sstart  send ...
            A = cols[0]
            B = cols[1]

            sA = sample_from_contig(A)
            sB = sample_from_contig(B)

            if (sA in EXCLUDE_SAMPLES) or (sB in EXCLUDE_SAMPLES):
                skipped += 1
                continue

            fout.write(line + "\n")
            kept += 1

    print(f"Kept lines: {kept}")
    print(f"Skipped lines: {skipped}")
    print(f"Wrote: {OUT_BLAST_GZ}")


if __name__ == "__main__":
    main()
