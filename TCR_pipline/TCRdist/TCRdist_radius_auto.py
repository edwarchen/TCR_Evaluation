#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import traceback
import pandas as pd
import numpy as np
from collections import deque

# 按照半径渐增的情况计算簇的情况，并且序列包含两条及以上才认为是簇，计算不同半径下所有簇的总信息（top3簇的频率，top10簇的频率等）

#########################################################
#  聚类函数：按 radius 合并
#  使用 BFS 合并距离 <= radius 的序列为 cluster
#########################################################
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

            try:
                neighbors = dist_mat.loc[cur][dist_mat.loc[cur] <= radius].index
            except Exception:
                neighbors = []

            for nb in neighbors:
                if nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
                    cluster.append(nb)

        clusters.append(cluster)

    return clusters

#########################################################
#  计算每个 radius 的 summary（统一簇定义：len >= 2）
#  并统计孤立序列数量 n_singletons
#########################################################
def summarize_clusters(clusters, freq_map):
    import numpy as np

    # -------------------------------
    # 1. 筛选有效簇（长度 >= 2）和孤立序列（长度=1）
    # -------------------------------
    valid_clusters = [cl for cl in clusters if len(cl) >= 2]
    singleton_clusters = [cl for cl in clusters if len(cl) == 1]  # 孤立序列
    n_singletons = len(singleton_clusters)

    # 若无有效簇 → 返回基本信息
    if len(valid_clusters) == 0:
        total_sequences = sum(freq_map.values())
        return {
            "n_sequence": int(total_sequences),        # 总序列数（包括孤立序列）
            "n_clusters": 0,                           # 有效簇数量（len>=2）
            "avg_cluster_size": 0.0,                   # 平均簇大小
            "largest_cluster_freq": 0.0,               # 最大簇频率（归一化）
            "largest_cluster_size": 0,                 # 最大簇大小
            "effective_cluster_number": 0.0,           # 有效簇数指标
            "gini_index": 0.0,                         # Gini 系数
            "top3_cluster_freq_sum": 0.0,              # top3簇频率和
            "top10_cluster_freq_sum": 0.0,             # top10簇频率和
            "n_singletons": n_singletons,              # 孤立序列数
        }

    # -------------------------------
    # 2. 计算有效簇的频率和大小
    # -------------------------------
    cluster_freqs = [sum(freq_map.get(x, 0.0) for x in cl) for cl in valid_clusters]
    cluster_sizes = [len(cl) for cl in valid_clusters]

    total_sequences = sum([len(cl) for cl in valid_clusters]) + n_singletons
    n_clusters = len(valid_clusters)
    avg_size = np.mean(cluster_sizes) if cluster_sizes else 0.0

    # 总频率用于归一化
    total_freq = sum(freq_map.values())
    if total_freq == 0:
        total_freq = 1.0

    # -------------------------------
    # 3. 计算 Gini 系数（衡量簇频率分布不均衡程度）
    # -------------------------------
    arr = np.array(cluster_freqs, dtype=float)
    if arr.sum() == 0:
        gini = 0.0
    else:
        arr_sorted = np.sort(arr)
        n = len(arr_sorted)
        gini = (
            2.0 * np.sum((np.arange(1, n + 1)) * arr_sorted)
            / (n * arr_sorted.sum())
            - (n + 1) / n
        )

    # -------------------------------
    # 4. top3 / top10 簇频率和（归一化）
    # -------------------------------
    valid_arr_desc = np.sort(arr)[::-1]  # 降序排列
    top3 = float(valid_arr_desc[:3].sum()) / total_freq
    top10 = float(valid_arr_desc[:10].sum()) / total_freq

    # 最大簇频率（归一化）和大小
    largest_cluster_freq_norm = float(max(valid_arr_desc)) / total_freq
    largest_cluster_size = int(max(cluster_sizes))

    # -------------------------------
    # 5. Effective cluster number（衡量簇均匀度）
    # -------------------------------
    sum_sq = sum(f * f for f in cluster_freqs)
    eff_clust = 1.0 / sum_sq if sum_sq > 0 else float("inf")

    # -------------------------------
    # 6. 返回结果
    # -------------------------------
    return {
        "n_sequence": int(total_sequences),                # 总序列数（包括孤立）
        "n_clusters": int(n_clusters),                     # 有效簇数量（len>=2）
        "avg_cluster_size": float(avg_size),              # 平均簇大小
        "largest_cluster_freq": float(largest_cluster_freq_norm), # 最大簇频率归一化
        "largest_cluster_size": int(largest_cluster_size),        # 最大簇大小
        "effective_cluster_number": float(eff_clust),     # 有效簇数指标
        "gini_index": float(gini),                        # Gini 系数
        "top3_cluster_freq_sum": float(top3),             # top3簇频率和
        "top10_cluster_freq_sum": float(top10),           # top10簇频率和
        "n_singletons": n_singletons,                     # 孤立序列数
    }

#########################################################
# 自动递增 radius，直到所有有效簇合并成一个
#########################################################
def auto_radius_until_single(dist_mat, freq_map, out_csv):
    results = []
    radius = 1

    while True:
        clusters = cluster_by_radius(dist_mat, radius)
        summary = summarize_clusters(clusters, freq_map)
        summary["radius"] = radius

        print(
            f"[INFO] radius={radius}, "
            f"n_clusters={summary['n_clusters']}, "
            f"n_singletons={summary['n_singletons']}"
        )
        results.append(summary)

        # ✅ 正确的停止条件：
        # 1 个有效簇 + 没有孤立序列
        if summary["n_clusters"] == 1 and summary["n_singletons"] == 0:
            print(
                f"[OK] radius={radius} 时，"
                f"所有序列进入同一个簇（无孤立序列），停止递增"
            )
            break

        radius += 1

    df = pd.DataFrame(results)
    df = df[
        [
            "radius",
            "n_sequence",
            "n_clusters",
            "avg_cluster_size",
            "largest_cluster_freq",
            "largest_cluster_size",
            "effective_cluster_number",
            "gini_index",
            "top3_cluster_freq_sum",
            "top10_cluster_freq_sum",
            "n_singletons",
        ]
    ]

    df.to_csv(out_csv, index=False)
    print(f"[OK] 写出 {out_csv}")


#########################################################
# 批处理流程
#########################################################
def process_batch(batch_dir):
    dist_files = []
    for root, _, files in os.walk(batch_dir):
        for f in files:
            if f.endswith("_pairwise_tcrdist.csv"):
                dist_files.append(os.path.join(root, f))

    if not dist_files:
        print(f"[WARN] 批次目录 {batch_dir} 下没有找到 *_pairwise_tcrdist.csv 文件")
        return

    dist_files.sort()

    for pairwise_file in dist_files:
        try:
            basename = os.path.basename(pairwise_file)
            sample_name = basename.replace("_pairwise_tcrdist.csv", "")

            print(f"\n===== 处理样本：{sample_name} =====")

            # 对应的 idmap
            idmap_file = pairwise_file.replace(
                "_pairwise_tcrdist.csv", "_unique_id_map.csv"
            )
            if not os.path.exists(idmap_file):
                print(f"[ERROR] 缺少 unique_id_map 文件：{idmap_file}")
                continue

            # 读取 pairwise matrix
            dist_mat = pd.read_csv(pairwise_file, index_col=0)

            # 读取 freq map
            df_idmap = pd.read_csv(idmap_file)
            if "unique_id" not in df_idmap.columns or "freq" not in df_idmap.columns:
                print("[ERROR] unique_id_map.csv 缺少必要字段 unique_id / freq")
                continue

            freq_map = dict(zip(df_idmap["unique_id"], df_idmap["freq"]))

            # 输出文件名
            out_csv = os.path.join(
                os.path.dirname(pairwise_file),
                f"{sample_name}_radius_summary.csv",
            )

            auto_radius_until_single(dist_mat, freq_map, out_csv)

        except Exception:
            traceback.print_exc()
            continue

#########################################################
# 主程序入口
#########################################################
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法：python TCRdist_radius_auto.py <batch_path>")
        sys.exit(1)

    batch_dir = sys.argv[1]
    process_batch(batch_dir)
print('[Done]')
