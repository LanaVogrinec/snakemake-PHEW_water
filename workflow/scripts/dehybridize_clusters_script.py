#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from collections import defaultdict

# Realms you want to recognize explicitly (optional; unknowns handled anyway)
REAL_REALMS = {
    "Duplodnaviria",
    "Monodnaviria",
    "Varidnaviria",
    "Riboviria",
    "Ribozyviria",
    "Adnaviria",
}

UNKNOWN_KEYWORDS = ["unclassified", "environmental", "uncultured", "sample"]


def is_unknown(val: str) -> bool:
    if val is None:
        return True
    v = str(val).strip().lower()
    if not v:
        return True
    return any(k in v for k in UNKNOWN_KEYWORDS)


def parse_lineage(lineage: str):
    """
    Extract realm/kingdom/phylum/class/order/family/genus from an ICTV-like lineage string.
    Returns dict with keys: realm, kingdom, phylum, class, order, family, genus
    """
    out = {"realm": None, "kingdom": None, "phylum": None, "class": None, "order": None, "family": None, "genus": None}
    if lineage is None or pd.isna(lineage):
        return out

    parts = [p.strip() for p in str(lineage).split(";") if p.strip()]

    # realm: first "-_" entry that is not Viruses
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
    """
    Normalize rank values so missing/unclassified collapse to '<rank>_unknown'.
    For realm, optionally restrict to REAL_REALMS (else unknown).
    """
    if val is None or pd.isna(val) or is_unknown(val):
        return f"{rank_name}_unknown"

    v = str(val).strip()
    if rank_name == "realm":
        # keep only recognized realms, everything else goes to unknown bucket
        if v in REAL_REALMS:
            return v
        return "realm_unknown"
    return v


def cluster_is_hybrid(contigs, contig_ranks):
    """
    Hybrid category simplified:
    0 = all contigs share the same value at each of realm..order (unknown bucket allowed)
    1 = any of realm..order differs
    Only meaningful for non-singletons.
    """
    ranks_to_check = ["realm", "kingdom", "phylum", "class", "order"]

    for r in ranks_to_check:
        s = set(contig_ranks[c][r] for c in contigs)
        if len(s) > 1:
            return 1
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="clusters_with_taxonomy.tsv (must include contig_id, cluster_id, lineage, singleton).")
    ap.add_argument("--final_assign", required=True, help="Output: contig_id, final_cluster_id")
    ap.add_argument("--progress", required=True, help="Output progress summary TSV")
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    required = {"contig_id", "cluster_id", "lineage", "singleton"}
    if not required.issubset(df.columns):
        raise ValueError(f"Input must contain columns: {sorted(required)}. Found: {list(df.columns)}")

    # Build per-contig normalized taxonomy ranks (once)
    contig_ranks = {}
    for _, row in df[["contig_id", "lineage"]].drop_duplicates().iterrows():
        c = row["contig_id"]
        tax = parse_lineage(row["lineage"])
        contig_ranks[c] = {
            "realm":   norm_rank_value("realm", tax["realm"]),
            "kingdom": norm_rank_value("kingdom", tax["kingdom"]),
            "phylum":  norm_rank_value("phylum", tax["phylum"]),
            "class":   norm_rank_value("class", tax["class"]),
            "order":   norm_rank_value("order", tax["order"]),
        }

    # Current assignment starts as original clusters
    current_assignment = dict(zip(df["contig_id"], df["cluster_id"]))

    # Build current clusters mapping: cluster -> list(contigs)
    clusters = defaultdict(list)
    for contig, cl in current_assignment.items():
        clusters[cl].append(contig)

    # Identify which original clusters are non-singletons (from file, fastest & consistent)
    # singleton column is per row; take cluster-level min/max
    cluster_singleton_flag = df.groupby("cluster_id")["singleton"].max().to_dict()

    def count_hybrids(current_clusters):
        n = 0
        for cl, contigs in current_clusters.items():
            if cluster_singleton_flag.get(cl, 0) == 1:
                continue
            if len(contigs) <= 1:
                continue
            if cluster_is_hybrid(contigs, contig_ranks) == 1:
                n += 1
        return n

    # Progress tracking:
    # rows: total_clusters, realm, kingdom, phylum, class, order
    progress_rows = []

    initial_total = len(clusters)
    initial_hybrid = count_hybrids(clusters)
    progress_rows.append(("total_clusters", initial_total, initial_hybrid))

    # Splitting sequence
    steps = [
        ("realm", "r"),
        ("kingdom", "k"),
        ("phylum", "p"),
        ("class", "c"),
        ("order", "o"),
    ]

    for rank_name, suffix in steps:
        # Find hybrid clusters among non-singletons
        to_split = {}
        for cl, contigs in clusters.items():
            if cluster_singleton_flag.get(cl, 0) == 1:
                continue
            if len(contigs) <= 1:
                continue
            if cluster_is_hybrid(contigs, contig_ranks) == 1:
                to_split[cl] = contigs

        # If nothing to split, just record progress and continue
        if not to_split:
            progress_rows.append((rank_name, len(clusters), 0))
            continue

        # Perform splitting for those clusters; keep others unchanged
        new_clusters = {}
        # first carry over clusters that are not being split
        for cl, contigs in clusters.items():
            if cl not in to_split:
                new_clusters[cl] = contigs

        # split each hybrid cluster by the given rank
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

        # update clusters mapping
        clusters = new_clusters

        # update singleton flags for newly created clusters
        # (new ones cannot be singletons-only by definition, but may become size 1 after splitting)
        # We compute fresh from clusters sizes for correctness.
        cluster_singleton_flag = {cl: (1 if len(contigs) == 1 else 0) for cl, contigs in clusters.items()}

        # count remaining hybrid clusters after this step
        remaining_hybrid = count_hybrids(clusters)
        progress_rows.append((rank_name, len(clusters), remaining_hybrid))

    # Write final assignments (includes singletons, unchanged or not)
    out_assign = pd.DataFrame({
        "contig_id": list(current_assignment.keys()),
        "cluster_id": list(current_assignment.values())
    }).sort_values(["cluster_id", "contig_id"])
    out_assign.to_csv(args.final_assign, sep="\t", index=False)

    # Write progress summary in your requested format:
    # step | total_clusters_after_step | hybrid_clusters_after_step
    prog_df = pd.DataFrame(progress_rows, columns=["step", "total_clusters", "hybrid_clusters"])
    prog_df.to_csv(args.progress, sep="\t", index=False)


if __name__ == "__main__":
    main()
