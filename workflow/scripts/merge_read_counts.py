#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--viral-reps", required=True)
    p.add_argument("--mode", required=True)

    p.add_argument("--sample-list", nargs="*")
    p.add_argument("--nki-files", nargs="*")
    p.add_argument("--carry-files", nargs="*")

    p.add_argument("--out")
    p.add_argument("--out-nki")
    p.add_argument("--out-carry")

    return p.parse_args()


def clean_id(x):
    return str(x).strip().replace("\r", "").replace("\t", "").replace(" ", "")


def load_coverm(file):
    df = pd.read_csv(file, sep="\t")
    df.columns = df.columns.str.strip()

    read_cols = [c for c in df.columns if "Read" in c and "Count" in c]
    if not read_cols:
        raise ValueError(f"No Read Count column in file: {file}")

    read_col = read_cols[0]

    sample = file.split("/")[-1].split("_coverm")[0].replace("08_", "")

    out = df[["Contig", read_col]].copy()
    out.columns = ["rep_contig", sample]

    # IMPORTANT: only clean whitespace, DO NOT modify biological IDs
    out["rep_contig"] = out["rep_contig"].apply(clean_id)

    return out.set_index("rep_contig")


def merge(files):
    dfs = [load_coverm(f) for f in files]

    out = dfs[0]
    for df in dfs[1:]:
        out = out.join(df, how="outer")

    return out


def main():
    args = parse_args()

    reps = pd.read_csv(args.viral_reps, sep="\t")
    reps["rep_contig"] = reps["rep_contig"].apply(clean_id)

    base = reps.set_index("rep_contig")

    if args.mode == "samples":
        merged = merge(args.sample_list)
        out = base.join(merged, how="left").fillna(0)

        out = out.reset_index()
        out.to_csv(args.out, sep="\t", index=False)

    elif args.mode == "controls":

        if args.nki_files:
            nki = merge(args.nki_files)
            nki_out = base.join(nki, how="left").fillna(0)
            nki_out = nki_out.reset_index()
            nki_out.to_csv(args.out_nki, sep="\t", index=False)

        if args.carry_files:
            carry = merge(args.carry_files)
            carry_out = base.join(carry, how="left").fillna(0)
            carry_out = carry_out.reset_index()
            carry_out.to_csv(args.out_carry, sep="\t", index=False)

    else:
        raise ValueError(f"Unknown mode: {args.mode}")


if __name__ == "__main__":
    main()