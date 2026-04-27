#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd
from collections import deque
from itertools import combinations

EPS = 1e-8

# =========================
# 固定的方差（来自你的随机抽样结果）
# =========================
VAR_DICT = {
    "TTT": 0.036424868310354075,
    "DDD": 0.0065059749829047475,
    "AAA": 0.0005595186433243391
}

ROOT_DIR = "/haplox/users/xuliu/TCR_Project/"

# =========================
# 聚类：基于半径的连通分量
# =========================
def cluster_by_radius(dist_mat, radius):
    ids = dist_mat.index.tolist()
    visited = set()
    clusters = []
    cid = 0

    for seq in ids:
        if seq in visited:
            continue
        queue = deque([seq])
        visited.add(seq)
        members = [seq]

        while queue:
            cur = queue.popleft()
            neighbors = dist_mat.loc[cur][dist_mat.loc[cur] <= radius].index
            for n in neighbors:
                if n not in visited:
                    visited.add(n)
                    members.append(n)
                    queue.append(n)

        clusters.append((cid, members))
        cid += 1

    rows = []
    for cid, seqs in clusters:
        for s in seqs:
            rows.append({"unique_id": s, "cluster_id": cid})

    return pd.DataFrame(rows)


# =========================
# TTT
# =========================
def compute_TTT(cluster_df, freq_df):
    df = cluster_df.merge(freq_df, on="unique_id", how="left")
    valid = df.groupby("cluster_id").filter(lambda x: len(x) >= 2)
    if valid.empty:
        return np.nan
    return valid["freq_norm"].sum()


# =========================
# DDD
# =========================
def compute_DDD(cluster_df, freq_df):
    df = cluster_df.merge(freq_df, on="unique_id", how="left")
    valid = df.groupby("cluster_id").filter(lambda x: len(x) >= 2)
    if valid.empty:
        return np.nan
    return valid.groupby("cluster_id")["freq_norm"].sum().max()


# =========================
# AAA
# =========================
def compute_AAA(cluster_df, dist_mat):
    valid = cluster_df.groupby("cluster_id").filter(lambda x: len(x) >= 2)
    clusters = valid["cluster_id"].unique()
    if len(clusters) <= 1:
        return np.nan

    cluster_seqs = {
        c: valid.loc[valid["cluster_id"] == c, "unique_id"].values
        for c in clusters
    }

    inter_dists = []
    for c1, c2 in combinations(clusters, 2):
        d = dist_mat.loc[
            cluster_seqs[c1],
            cluster_seqs[c2]
        ].values.mean()
        inter_dists.append(d)

    inter_dist = np.mean(inter_dists)
    max_dist = dist_mat.values.max()
    return 1.0 - inter_dist / max_dist


# =========================
# CI：方差比例加权（动态指标）
# =========================
def compute_CI(ttt, ddd, aaa):
    metrics = {}
    if not np.isnan(ttt):
        metrics["TTT"] = ttt
    if not np.isnan(ddd):
        metrics["DDD"] = ddd
    if not np.isnan(aaa):
        metrics["AAA"] = aaa

    if len(metrics) == 0:
        return np.nan, {}

    vars_used = {k: VAR_DICT[k] + EPS for k in metrics}
    total_var = sum(vars_used.values())
    weights = {k: v / total_var for k, v in vars_used.items()}

    ci = sum(weights[k] * metrics[k] for k in metrics)
    return ci, weights


# =========================
# 主流程
# =========================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", required=True)
    parser.add_argument("--radius", type=float, required=True)
    args = parser.parse_args()

    batch_dir = os.path.join(ROOT_DIR, args.batch)

    samples = [
        d for d in os.listdir(batch_dir)
        if os.path.isdir(os.path.join(batch_dir, d))
    ]

    records = []

    for sample in samples:
        out_dir = os.path.join(
            batch_dir, sample, "Map_Clone_Analysis", "TCRdist_output"
        )

        dist_file = os.path.join(
            out_dir, f"{sample}.clonotypes.TRB_pairwise_tcrdist.csv"
        )
        freq_file = os.path.join(
            out_dir, f"{sample}.clonotypes.TRB_unique_id_map.csv"
        )

        if not (os.path.exists(dist_file) and os.path.exists(freq_file)):
            continue

        dist_mat = pd.read_csv(dist_file, index_col=0)

        freq_df = pd.read_csv(freq_file)
        freq_df = freq_df[["unique_id", "freq"]]
        freq_df["freq_norm"] = freq_df["freq"] / freq_df["freq"].sum()

        cluster_df = cluster_by_radius(dist_mat, args.radius)

        ttt = compute_TTT(cluster_df, freq_df)
        ddd = compute_DDD(cluster_df, freq_df)
        aaa = compute_AAA(cluster_df, dist_mat)

        ci, weights = compute_CI(ttt, ddd, aaa)

        records.append({
            "Sample_ID": sample,
            "radius": args.radius,
            "TTT": ttt,
            "DDD": ddd,
            "AAA": aaa,
            "TCR_CI": ci,
            "weights": weights
        })

    out_df = pd.DataFrame(records)
    out_path = os.path.join(
        ROOT_DIR, args.batch, f"TCRdist/TCR_CI_radius_{args.radius}.csv"
    )
    out_df.to_csv(out_path, index=False)
    print(f"[DONE] Saved {out_path}")


if __name__ == "__main__":
    main()
