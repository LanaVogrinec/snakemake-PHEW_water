#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from collections import defaultdict
from datetime import datetime

REAL_REALMS = {
    "Duplodnaviria",
    "Monodnaviria",
    "Varidnaviria",
    "Riboviria",
    "Ribozyviria",
    "Adnaviria",
}

UNKNOWN_KEYWORDS = ["unclassified", "environmental", "uncultured", "sample"]

def parse_args():
    p = argparse.ArgumentParser(description="Split hybrid clusters based on taxonomy ranks and track progress")
    p.add_argument("--input", required=True, help="clusters_with_taxonomy.tsv (must include contig_id, cluster_id, lineage, singleton)")
    p.add_argument("--final_assign", required=True, help="Output: contig_id, final_cluster_id")
    p.add_argument("--progress", required=True, help="Output progress summary TSV")
    return p.parse_args()

def log(msg):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", flush=True)

def is_unknown(val: str) -> bool:
    if val is None:
        return True
    v = str(val).strip().lower()
    if not v:
        return True
    return any(k in v for k in UNKNOWN_KEYWORDS)

def parse_lineage(lineage: str):
    out = {"realm": None, "kingdom": None, "phylum": None, "class": None, "order": None, "family": None, "genus": None}
    if lineage is None or pd.isna(lineage):
        return out
    parts = [p.strip() for p in str(lineage).split(";") if p.strip()]
    for p in parts:
        if p.startswith("-_"):
            name = p.replace("-_", "").strip()
            if name != "Viruses":
                out["realm"] = name
                break
    for p in parts:
        if p.startswith("k_"):
            out["kingdom"] = p.replace("k_", "").strip()
        elif p.startswith("p_"):
            out["phylum"] = p.replace("p_", "").strip()
        elif p.startswith("c_"):
            out["class"] = p.replace("c_", "").strip()
        elif p.startswith("o_"):
            out["order"] = p.replace("o_", "").strip()
        elif p.startswith("f_"):
            out["family"] = p.replace("f_", "").strip()
        elif p.startswith("g_"):
            out["genus"] = p.replace("g_", "").strip()
    return out

def norm_rank_value(rank_name: str, val: str) -> str:
    if val is None or pd.isna(val) or is_unknown(val):
        return f"{rank_name}_unknown"
    v = str(val).strip()
    if rank_name == "realm":
        return v if v in REAL_REALMS else "realm_unknown"
    return v

def cluster_is_hybrid(contigs, contig_ranks):
    ranks_to_check = ["realm", "kingdom", "phylum", "class", "order"]
    for r in ranks_to_check:
        if len({contig_ranks[c][r] for c in contigs}) > 1:
            return 1
    return 0

def main():
    args = parse_args()
    log(f"Reading input table: {args.input}")
    df = pd.read_csv(args.input, sep="\t")

    required = {"contig_id", "cluster_id", "lineage", "singleton"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input must contain columns: {sorted(required)}")

    log("Normalizing contig taxonomy")
    contig_ranks = {}
    for _, row in df[["contig_id", "lineage"]].drop_duplicates().iterrows():
        c = row["contig_id"]
        tax = parse_lineage(row["lineage"])
        contig_ranks[c] = {r: norm_rank_value(r, tax[r]) for r in ["realm","kingdom","phylum","class","order"]}

    current_assignment = dict(zip(df["contig_id"], df["cluster_id"]))
    clusters = defaultdict(list)
    for contig, cl in current_assignment.items():
        clusters[cl].append(contig)

    cluster_singleton_flag = df.groupby("cluster_id")["singleton"].max().to_dict()

    def count_hybrids(current_clusters):
        n = 0
        for cl, contigs in current_clusters.items():
            if cluster_singleton_flag.get(cl, 0) == 1 or len(contigs) <= 1:
                continue
            if cluster_is_hybrid(contigs, contig_ranks):
                n += 1
        return n

    progress_rows = []
    progress_rows.append(("total_clusters", len(clusters), count_hybrids(clusters)))

    steps = [("realm", "r"), ("kingdom", "k"), ("phylum", "p"), ("class", "c"), ("order", "o")]

    for rank_name, suffix in steps:
        to_split = {cl: contigs for cl, contigs in clusters.items() 
                    if cluster_singleton_flag.get(cl, 0) == 0 and len(contigs) > 1 and cluster_is_hybrid(contigs, contig_ranks)}
        if not to_split:
            progress_rows.append((rank_name, len(clusters), 0))
            continue

        new_clusters = {cl: contigs for cl, contigs in clusters.items() if cl not in to_split}
        for old_cl, contigs in to_split.items():
            groups = defaultdict(list)
            for c in contigs:
                groups[contig_ranks[c][rank_name]].append(c)
            sub_idx = 1
            for _, group_contigs in groups.items():
                new_cl = f"{old_cl}_{suffix}{sub_idx}"
                new_clusters[new_cl] = group_contigs
                for c in group_contigs:
                    current_assignment[c] = new_cl
                sub_idx += 1

        clusters = new_clusters
        cluster_singleton_flag = {cl: (1 if len(contigs) == 1 else 0) for cl, contigs in clusters.items()}
        progress_rows.append((rank_name, len(clusters), count_hybrids(clusters)))

    log(f"Writing final cluster assignments: {args.final_assign}")
    out_assign = pd.DataFrame({
        "contig_id": list(current_assignment.keys()),
        "cluster_id": list(current_assignment.values())
    }).sort_values(["cluster_id", "contig_id"])
    out_assign.to_csv(args.final_assign, sep="\t", index=False)

    log(f"Writing progress summary: {args.progress}")
    prog_df = pd.DataFrame(progress_rows, columns=["step", "total_clusters", "hybrid_clusters"])
    prog_df.to_csv(args.progress, sep="\t", index=False)

    log("Cluster splitting completed")

if __name__ == "__main__":
    main()
