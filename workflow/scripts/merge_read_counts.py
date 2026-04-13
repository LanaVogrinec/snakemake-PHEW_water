#!/usr/bin/env python3

import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--viral-reps", required=True)
    p.add_argument("--taxonomy", required=True)

    p.add_argument("--mode", required=True)

    p.add_argument("--sample-list", nargs="*")
    p.add_argument("--nki-files", nargs="*")
    p.add_argument("--carry-files", nargs="*")

    p.add_argument("--out")
    p.add_argument("--out-nki")
    p.add_argument("--out-carry")

    return p.parse_args()


def load_one(file):
    df = pd.read_csv(file, sep="\t")

    sample = file.split("/")[-1].split("_coverm")[0].replace("08_", "")
    col = [c for c in df.columns if "Read Count" in c][0]

    out = df[["Contig", col]].copy()
    out.columns = ["Contig", sample]
    out["Contig"] = out["Contig"].str.replace(r"_cluster_.*$", "", regex=True)

    return out.set_index("Contig")


def merge(files):
    if not files:
        return None

    dfs = [load_one(f) for f in files]

    out = dfs[0]
    for df in dfs[1:]:
        out = out.join(df, how="outer")

    return out


def main():
    args = parse_args()

    reps = pd.read_csv(args.viral_reps, sep="\t")
    tax = pd.read_csv(args.taxonomy, sep="\t", header=None)
    tax.columns = ["rep_contig", "taxonomy"]

    base = reps.merge(tax, on="rep_contig", how="left")

    # -------------------
    # SAMPLES MODE
    # -------------------
    if args.mode == "samples" and args.sample_list:
        merged = merge(args.sample_list)
        out = base.join(merged, on="rep_contig")
        out.to_csv(args.out, sep="\t", index=False)

    # -------------------
    # CONTROLS MODE
    # -------------------
    if args.mode == "controls":

        if args.nki_files:
            nki = merge(args.nki_files)
            out_nki = base.join(nki, on="rep_contig")
            out_nki.to_csv(args.out_nki, sep="\t", index=False)

        if args.carry_files:
            carry = merge(args.carry_files)
            out_carry = base.join(carry, on="rep_contig")
            out_carry.to_csv(args.out_carry, sep="\t", index=False)


if __name__ == "__main__":
    main()