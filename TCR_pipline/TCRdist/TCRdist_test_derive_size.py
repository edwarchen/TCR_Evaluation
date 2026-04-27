#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Phase1: 派生样本稳定性扫描脚本（多进程优化）
目标：
1. 探索不同派生样本数量(N_derive)下，TTT、DDD、AAA的稳定性
2. 输出指标均值、标准差和CV，用于确定稳定的派生样本数量

输出文件：
1. derived_raw_metrics.tsv：每个派生样本的指标
2. derived_stability_summary.tsv：每个N_derive下指标均值、标准差、CV
"""

import os
import pandas as pd
import numpy as np
from collections import deque
from itertools import combinations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

EPS = 1e-8

# -------------------------------
# 聚类函数：基于半径的连通分量
# -------------------------------
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

# -------------------------------
# 指标计算
# -------------------------------
def compute_metrics(cluster_df, dist_mat, freq_df, aaa_const=None):
    """
    cluster_df: 当前簇分配
    dist_mat: 当前使用的距离矩阵
    freq_df: unique_id与freq_norm
    aaa_const: 原始样本AAA，不为空时派生样本直接使用
    """
    valid = cluster_df.groupby("cluster_id").filter(lambda x: len(x) >= 2)
    n_cluster = valid["cluster_id"].nunique()

    # Case1: 无有效簇
    if n_cluster == 0:
        return n_cluster, np.nan, np.nan, np.nan

    # 合并频率
    df = valid.merge(freq_df, on="unique_id", how="left")
    cluster_freq = df.groupby("cluster_id")["freq_norm"].sum()

    ttt = cluster_freq.sum()
    ddd = cluster_freq.max()

    # Case2: 仅1个簇
    if n_cluster == 1:
        aaa = aaa_const  # 如果派生样本传入原样本AAA
        return n_cluster, aaa, ttt, ddd

    # Case3: 多个簇 → AAA计算
    if aaa_const is not None:
        aaa = aaa_const
    else:
        cluster_seqs = {c: valid.loc[valid["cluster_id"]==c,"unique_id"].values
                        for c in valid["cluster_id"].unique()}
        inter_dists = []
        for c1, c2 in combinations(cluster_seqs.keys(), 2):
            inter_dists.append(dist_mat.loc[cluster_seqs[c1], cluster_seqs[c2]].values.mean())
        inter_dist = np.mean(inter_dists)
        max_dist = dist_mat.values.max()
        aaa = 1.0 - inter_dist / (max_dist + EPS)
    return n_cluster, aaa, ttt, ddd

# -------------------------------
# 样本处理函数（每个样本独立计算）
# -------------------------------
def process_sample(row, radius, derive_sizes):
    batch = row["batch_id"]
    sample = row["sample_id"]
    root = "/haplox/users/xuliu/TCR_Project/Results"
    out_dir = os.path.join(root, batch, sample, "Map_Clone_Analysis", "TCRdist_output")
    pre_file = os.path.join(out_dir, f"{sample}.clonotypes.TRB_preprocessed.tsv")
    dist_file = os.path.join(out_dir, f"{sample}.clonotypes.TRB_pairwise_tcrdist.csv")

    if not os.path.exists(pre_file) or not os.path.exists(dist_file):
        print(f"[WARN] Missing input for {sample}")
        return []

    pre = pd.read_csv(pre_file, sep="\t")
    dist_full = pd.read_csv(dist_file, index_col=0)

    pre["freq_norm"] = pre["count"] / pre["count"].sum()

    # 先计算原始样本AAA
    cluster_orig = cluster_by_radius(dist_full, radius)
    _, aaa_const, _, _ = compute_metrics(cluster_orig, dist_full, pre)

    sample_records = []

    for N in derive_sizes:
        for i in range(N):
            derive = pre.sample(frac=1.0, replace=True)
            used_ids = derive["unique_id"].values
            dist_mat = dist_full.loc[used_ids, used_ids]
            freq_df = derive[["unique_id","freq_norm"]]
            cluster_df = cluster_by_radius(dist_mat, radius)
            n_cluster, aaa, ttt, ddd = compute_metrics(cluster_df, dist_mat, freq_df, aaa_const)

            sample_records.append({
                "batch": batch,
                "sample": sample,
                "N_derive": N,
                "derive_id": i,
                "n_cluster": n_cluster,
                "AAA": aaa,
                "TTT": ttt,
                "DDD": ddd
            })
    print(f"[DONE] {sample}")
    return sample_records

# -------------------------------
# 主流程
# -------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--phase_samples", type=int, default=3,
                        help="Phase1中处理样本数量，默认3个")
    parser.add_argument("--max_workers", type=int, default=3,
                        help="同时运行样本进程数")
    args = parser.parse_args()

    sample_list = "/haplox/users/xuliu/TCR_Project/SampleList/Derived_analysis.tsv"
    samples = pd.read_csv(sample_list, sep="\t")
    phase_samples = samples.head(args.phase_samples)

    derive_sizes = [5,10,20,30,50,80,100]
    all_records = []

    # 多进程处理样本
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {executor.submit(process_sample, row, args.radius, derive_sizes): row["sample_id"]
                   for _, row in phase_samples.iterrows()}
        for fut in as_completed(futures):
            records = fut.result()
            all_records.extend(records)

    # 保存原始派生指标
    raw_df = pd.DataFrame(all_records)
    raw_df.to_csv("derived_raw_metrics.tsv", sep="\t", index=False)

    # 汇总稳定性
    summary = (
        raw_df
        .melt(id_vars=["batch","sample","N_derive"], value_vars=["AAA","TTT","DDD"],
              var_name="metric", value_name="value")
        .groupby(["sample","N_derive","metric"])
        .agg(mean=("value","mean"), sd=("value","std"))
        .reset_index()
    )
    summary["cv"] = summary["sd"] / (summary["mean"].abs() + EPS)
    summary.to_csv("derived_stability_summary.tsv", sep="\t", index=False)

    print("[FINISHED] Phase1 Stability scan completed.")

if __name__ == "__main__":
    main()
