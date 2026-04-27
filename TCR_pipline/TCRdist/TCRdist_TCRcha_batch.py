#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
from collections import deque
from itertools import combinations

# -----------------------
# 多样性指标函数
# -----------------------
def shannon_index(freqs):
    freqs = np.array(freqs)
    freqs = freqs[freqs > 0]
    return -np.sum(freqs * np.log2(freqs)) if len(freqs) > 0 else 0.0

def simpson_index(freqs):
    freqs = np.array(freqs)
    return 1 - np.sum(freqs**2)

def evenness_index(freqs):
    n = len(freqs)
    return shannon_index(freqs)/np.log2(n) if n > 1 else 0.0

def clonality_index(freqs):
    return 1 - evenness_index(freqs)

def observed_diversity_mean(dist_mat, clusters):
    dists = []
    for cl in clusters:
        if len(cl) <= 1:
            continue  # 单序列 cluster 无法计算 pairwise distance
        pair_dists = [dist_mat.loc[i, j] for i, j in combinations(cl, 2)]
        dists.append(np.mean(pair_dists))
    if len(dists) == 0:
        return np.nan
    return np.mean(dists)

# -----------------------
# 聚类函数
# -----------------------
def cluster_by_radius(dist_mat, radius):
    ids = list(dist_mat.index)
    visited = set()
    clusters = []
    for seq in ids:
        if seq in visited:
            continue
        queue = deque([seq])
        cluster = [seq]
        visited.add(seq)
        while queue:
            cur = queue.popleft()
            neighbors = dist_mat.loc[cur][dist_mat.loc[cur] <= radius].index
            for nb in neighbors:
                if nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
                    cluster.append(nb)
        # 不再要求 len(cluster) >= 2
        clusters.append(cluster)
    return clusters

# -----------------------
# 主函数
# -----------------------
def main(batch_dir, radius):
    batch_tcrdist_dir = os.path.join(batch_dir, "TCRdist")
    os.makedirs(batch_tcrdist_dir, exist_ok=True)
    out_dir = os.path.join(batch_tcrdist_dir, f"cluster_metrics_R{radius}")
    os.makedirs(out_dir, exist_ok=True)

    metrics_all_samples = []

    for sample in os.listdir(batch_dir):
        sample_path = os.path.join(batch_dir, sample)

        if not os.path.isdir(sample_path) or sample == "TCRdist":
            continue

        tcr_output_dir = os.path.join(sample_path, "Map_Clone_Analysis", "TCRdist_output")
        if not os.path.isdir(tcr_output_dir):
            continue

        pairwise_file = None
        uid_map_file = None
        for f in os.listdir(tcr_output_dir):
            if "TRB_pairwise_tcrdist" in f and f.endswith(".csv"):
                pairwise_file = os.path.join(tcr_output_dir, f)
            elif "TRB_unique_id_map" in f and f.endswith(".csv"):
                uid_map_file = os.path.join(tcr_output_dir, f)
        if not pairwise_file or not uid_map_file:
            continue

        df_uid = pd.read_csv(uid_map_file)
        df_dist = pd.read_csv(pairwise_file, index_col=0)

        df_uid["freq_norm"] = df_uid["freq"] / df_uid["freq"].sum()
        clusters = cluster_by_radius(df_dist, radius)

        # 保存序列对应 cluster
        cluster_assign = []
        for cid, clust in enumerate(clusters, 1):
            for seq in clust:
                cluster_assign.append({
                    "unique_id": seq,
                    "cluster_id": f"C{cid}"
                })
        df_cluster_map = pd.DataFrame(cluster_assign)
        df_cluster_map = df_cluster_map.merge(df_uid, on="unique_id", how="left")

        sample_out_dir = os.path.join(out_dir, sample)
        os.makedirs(sample_out_dir, exist_ok=True)
        seq_cluster_file = os.path.join(sample_out_dir, "sequence_cluster.tsv")
        df_cluster_map[["unique_id","cluster_id","freq_norm","v_b_gene","j_b_gene"]].to_csv(seq_cluster_file, sep="\t", index=False)

        # 计算 cluster 频率
        cluster_freqs = []
        for clust in clusters:
            freq_sum = df_uid[df_uid["unique_id"].isin(clust)]["freq"].sum()
            cluster_freqs.append(freq_sum / df_uid["freq"].sum())

        metrics = {
            "sample": sample,
            "n_clusters": len(clusters),
            "shannon": shannon_index(cluster_freqs),
            "simpson": simpson_index(cluster_freqs),
            "evenness": evenness_index(cluster_freqs),
            "clonality": clonality_index(cluster_freqs),
            "InternalDistance_mean": observed_diversity_mean(df_dist, clusters)
        }
        metrics_all_samples.append(metrics)

    metrics_df = pd.DataFrame(metrics_all_samples)
    metrics_file = os.path.join(out_dir, "sample_level_cluster_metrics.tsv")
    metrics_df.to_csv(metrics_file, sep="\t", index=False)
    print(f"OUTPUT: {metrics_file}")

# -----------------------
# CLI
# -----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute cluster-level TCR metrics for batch")
    parser.add_argument("--batch_dir", required=True, help="Batch directory containing sample folders")
    parser.add_argument("--radius", required=True, type=float, help="Radius for clustering")
    args = parser.parse_args()
    main(args.batch_dir, args.radius)
