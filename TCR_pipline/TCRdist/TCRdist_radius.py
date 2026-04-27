#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict, deque

# python TCRdist_radius.py \
#  convert.sample.tcrdist_matrix.csv \
#  convert.sample_unique_id_map.csv \
#  1,3,5,7,9,10,15,20,25,30 \
#  /haplox/users/xuliu/TCR_Project/Results/250716_TCRcell_supple/S055_SZ20250313034TUT-8_ttdna_genome_235590/Map_Clone_Analysis/TCRdist_output/


############################################
#  基础函数：给定 radius，做按距离的聚类（相似度聚合）
############################################
def cluster_by_radius(dist_mat, radius):
    """
    dist_mat: pandas DataFrame，行为 unique_id，列为 unique_id
    radius: 半径阈值
    return: list of clusters，每个 cluster 是一个序列 ID 列表
    """

    ids = list(dist_mat.index)
    visited = set()
    clusters = []

    for seq in ids:
        if seq in visited:
            continue

        # BFS 聚类
        queue = deque([seq])
        cluster = [seq]
        visited.add(seq)

        while queue:
            cur = queue.popleft()
            # 找所有距离 <= radius 的邻居
            neighbors = dist_mat.loc[cur][dist_mat.loc[cur] <= radius].index

            for nb in neighbors:
                if nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
                    cluster.append(nb)

        clusters.append(cluster)

    return clusters


############################################
#  统计每个 radius 的 summary 信息
############################################
def summarize_clusters(clusters, freq_map):
    """
    clusters: list of cluster(list)
    freq_map: dict {unique_id: freq}
    """

    # 簇频率总和
    cluster_freqs = []
    cluster_sizes = []

    for cl in clusters:
        fs = sum(freq_map.get(x, 0) for x in cl)
        cluster_freqs.append(fs)
        cluster_sizes.append(len(cl))

    n_clusters = len(clusters)
    total_sequences = sum(cluster_sizes)

    # 最大簇
    largest_cluster_freq = max(cluster_freqs)
    largest_cluster_size = max(cluster_sizes)

    # 平均簇大小
    avg_size = total_sequences / n_clusters if n_clusters > 0 else 0

    # 有效簇数（1 / sum(freq^2)）
    eff_clust = 1 / sum(f * f for f in cluster_freqs)

    # Gini 系数（衡量聚类集中度）
    arr = np.array(cluster_freqs)
    arr_sorted = np.sort(arr)
    n = len(arr_sorted)
    gini = (2 * np.sum((np.arange(1, n + 1)) * arr_sorted) / (n * arr_sorted.sum())) - (n + 1) / n

    # Top3 和 Top10 频率
    arr_desc = arr_sorted[::-1]
    top3 = arr_desc[:3].sum()
    top10 = arr_desc[:10].sum()

    return {
        "n_sequence": total_sequences,
        "n_clusters": n_clusters,
        "avg_cluster_size": avg_size,
        "largest_cluster_freq": largest_cluster_freq,
        "largest_cluster_size": largest_cluster_size,
        "effective_cluster_number": eff_clust,
        "gini_index": gini,
        "top3_cluster_freq_sum": top3,
        "top10_cluster_freq_sum": top10,
    }


############################################
#  主函数：批量样本 or 单一样本
############################################
def main():
    if len(sys.argv) < 5:
        print("Usage:")
        print("  python TCRdist_radius.py <dist_matrix.csv> <id_map.csv> <radius_list> <output_dir>")
        print("  radius_list 例如: 1,3,5,7,10")
        sys.exit(1)

    dist_file = sys.argv[1]
    idmap_file = sys.argv[2]
    radius_list = [int(x) for x in sys.argv[3].split(",")]
    output_dir = sys.argv[4]

    os.makedirs(output_dir, exist_ok=True)

    print(f"读取距离矩阵: {dist_file}")
    dist_mat = pd.read_csv(dist_file, index_col=0)

    print(f"读取 ID map: {idmap_file}")
    idmap = pd.read_csv(idmap_file)

    if "unique_id" not in idmap.columns or "freq" not in idmap.columns:
        print("id_map.csv 必须包含 unique_id 和 freq 列")
        sys.exit(1)

    # ID 对齐
    seq_ids = dist_mat.index
    idmap = idmap.set_index("unique_id").loc[seq_ids]
    freq_map = idmap["freq"].to_dict()

    # 输出 summary
    summary_records = []

    for R in radius_list:
        print(f"计算 radius = {R} ...")

        clusters = cluster_by_radius(dist_mat, R)
        summary = summarize_clusters(clusters, freq_map)
        summary["radius"] = R

        summary_records.append(summary)

    out_csv = os.path.join(output_dir, "radius_summary.csv")
    pd.DataFrame(summary_records).to_csv(out_csv, index=False)

    print(f"\n已生成：{out_csv}")


if __name__ == "__main__":
    main()
