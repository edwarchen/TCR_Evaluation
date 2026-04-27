#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import numpy as np
from collections import deque
from itertools import combinations
import glob

#对批次下的样本，计算逐步增加半径后的TCR特征
# python TCR_Project/scripts/TCRdist/TCRdist_TCRcha_batch_multiradius.py --batch_dir /haplox/users/xuliu/TCR_Project/Results/250817_TCRcell/ --max_radius 300（最大径设置）

# =========================
# 聚类：半径法
# =========================
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

        clusters.append(cluster)

    return clusters


# =========================
# 多样性指标
# =========================
def shannon(freqs):
    freqs = np.array(freqs)
    freqs = freqs[freqs > 0]
    return -np.sum(freqs * np.log(freqs))

def simpson(freqs):
    freqs = np.array(freqs)
    return np.sum(freqs ** 2)

def evenness(shannon_val, n):
    if n <= 1:
        return np.nan
    return shannon_val / np.log(n)

def clonality(evenness_val):
    return 1 - evenness_val


def observed_diversity_mean(dist_mat, clusters):
    dists = []
    for cl in clusters:
        if len(cl) <= 1:
            continue
        pair_dists = [
            dist_mat.loc[i, j]
            for i, j in combinations(cl, 2)
        ]
        dists.append(np.mean(pair_dists))

    if len(dists) == 0:
        return np.nan
    return np.mean(dists)


# =========================
# 单样本处理
# =========================
def process_sample(sample_name, map_clone_dir, out_base, max_radius):

    tcrdist_out = os.path.join(map_clone_dir, "TCRdist_output")

    dist_files = glob.glob(
        os.path.join(tcrdist_out, "*.clonotypes.TRB_pairwise_tcrdist.csv")
    )
    freq_files = glob.glob(
        os.path.join(tcrdist_out, "*.clonotypes.TRB_unique_id_map.csv")
    )

    if len(dist_files) == 0 or len(freq_files) == 0:
        print(f"[跳过] {sample_name}：未找到 TCRdist 输入文件")
        return

    dist_file = dist_files[0]
    freq_file = freq_files[0]

    print(f"[处理] {sample_name}")
    print(f"  dist: {os.path.basename(dist_file)}")
    print(f"  freq: {os.path.basename(freq_file)}")

    dist_mat = pd.read_csv(dist_file, index_col=0)
    freq_df = pd.read_csv(freq_file)

    freq_map = dict(zip(freq_df["unique_id"], freq_df["freq"]))

    results = []

    for r in range(1, max_radius + 1):
        clusters = cluster_by_radius(dist_mat, r)

        stop = (len(clusters) == 1)

        cluster_freqs = [
            sum(freq_map.get(seq, 0) for seq in cl)
            for cl in clusters
        ]

        total = sum(cluster_freqs)
        if total == 0:
            continue

        cluster_freqs = [f / total for f in cluster_freqs]

        sh = shannon(cluster_freqs)
        sp = simpson(cluster_freqs)
        ev = evenness(sh, len(cluster_freqs))
        cl = clonality(ev)
        obs = observed_diversity_mean(dist_mat, clusters)

        results.append({
            "radius": r,
            "n_clusters": len(clusters),
            "shannon": sh,
            "simpson": sp,
            "evenness": ev,
            "clonality": cl,
            "observedDiversity_mean": obs
        })

        if stop:
            break

    out_dir = os.path.join(out_base, sample_name)
    os.makedirs(out_dir, exist_ok=True)

    out_file = os.path.join(out_dir, "metrics_by_radius.tsv")
    pd.DataFrame(results).to_csv(out_file, sep="\t", index=False)

    print(f"  输出 → {out_file}")


# =========================
# batch 主程序
# =========================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch_dir", required=True)
    parser.add_argument("--max_radius", type=int, default=300)
    args = parser.parse_args()

    batch_dir = args.batch_dir
    out_base = os.path.join(batch_dir, "TCRdist")
    os.makedirs(out_base, exist_ok=True)

    for d in os.listdir(batch_dir):
        sample_path = os.path.join(batch_dir, d)
        map_clone = os.path.join(sample_path, "Map_Clone_Analysis")

        if os.path.isdir(map_clone):
            process_sample(d, map_clone, out_base, args.max_radius)


if __name__ == "__main__":
    main()
