#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import traceback
import pandas as pd
import numpy as np
from collections import deque


# 按照给定的半径计算簇的情况，并且序列包含两条及以上才认为是簇，计算不同半径下所有簇的总信息（top3簇的频率，top10簇的频率等）
# usage:
# python TCR_Project/scripts/TCRdist/TCRdist_radius_batch.py /haplox/users/xuliu/TCR_Project/Results/250712_TCRcell/ 1,3,5,7,10

############################################
#  基础函数：给定 radius，做按距离的聚类（相似度聚合）
############################################
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
            # 找到距离 <= radius 的邻居（包含自身）
            # 使用 .values <= radius 速度稍好，但保持可读性用布尔索引
            try:
                neighbors = dist_mat.loc[cur][dist_mat.loc[cur] <= radius].index
            except Exception:
                # 如果行访问失败（极少数），跳过该节点
                neighbors = []
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
    cluster_freqs = []
    cluster_sizes = []

    for cl in clusters:
        fs = sum(freq_map.get(x, 0.0) for x in cl)
        cluster_freqs.append(float(fs))
        cluster_sizes.append(len(cl))

    if len(cluster_freqs) == 0:
        return {
            "n_sequence": 0,
            "n_clusters": 0,
            "avg_cluster_size": 0,
            "largest_cluster_freq": 0.0,
            "largest_cluster_size": 0,
            "effective_cluster_number": 0.0,
            "gini_index": 0.0,
            "top3_cluster_freq_sum": 0.0,
            "top10_cluster_freq_sum": 0.0,
        }

    n_clusters = len(clusters)
    total_sequences = sum(cluster_sizes)
    largest_cluster_freq = max(cluster_freqs) if cluster_freqs else 0.0
    largest_cluster_size = max(cluster_sizes) if cluster_sizes else 0
    avg_size = total_sequences / n_clusters if n_clusters > 0 else 0.0

    sum_sq = sum(f * f for f in cluster_freqs)
    eff_clust = 1.0 / sum_sq if sum_sq > 0 else float('inf')

    arr = np.array(cluster_freqs, dtype=float)
    total_freq = arr.sum() if arr.sum() > 0 else 1.0  # 防止除零

    if arr.sum() == 0:
        gini = 0.0
    else:
        arr_sorted = np.sort(arr)
        n = len(arr_sorted)
        gini = (2.0 * np.sum((np.arange(1, n + 1)) * arr_sorted) / (n * arr_sorted.sum())) - (n + 1) / n

    arr_desc = np.sort(arr)[::-1]
    # 归一化为占比
    largest_cluster_freq /= total_freq
    top3 = float(arr_desc[:3].sum()) / total_freq if arr_desc.size > 0 else 0.0
    top10 = float(arr_desc[:10].sum()) / total_freq if arr_desc.size > 0 else 0.0

    return {
        "n_sequence": int(total_sequences),
        "n_clusters": int(n_clusters),
        "avg_cluster_size": float(avg_size),
        "largest_cluster_freq": float(largest_cluster_freq),
        "largest_cluster_size": int(largest_cluster_size),
        "effective_cluster_number": float(eff_clust),
        "gini_index": float(gini),
        "top3_cluster_freq_sum": float(top3),
        "top10_cluster_freq_sum": float(top10),
    }


############################################
#  批量处理函数（对每个样本 try/except，确保处理全部样本）
############################################
def process_batch(batch_dir, radius_list):
    # 收集所有 _pairwise_tcrdist.csv 文件（包括子目录）
    dist_files = []
    for root, _, files in os.walk(batch_dir):
        for f in files:
            if f.endswith("_pairwise_tcrdist.csv"):
                dist_files.append(os.path.join(root, f))

    if not dist_files:
        print(f"[WARN] 在 {batch_dir} 没有找到任何 *_pairwise_tcrdist.csv 文件")
        return

    dist_files.sort()  # 稳定顺序（可选）

    for dist_file in dist_files:
        try:
            sample_name = os.path.basename(dist_file).replace("_pairwise_tcrdist.csv", "")
            idmap_file = dist_file.replace("_pairwise_tcrdist.csv", "_unique_id_map.csv")

            if not os.path.exists(idmap_file):
                print(f"[WARN] {sample_name} 缺少 id_map 文件（{idmap_file}），跳过该样本")
                continue

            print(f"\n[INFO] 处理样本: {sample_name}")
            out_dir = os.path.join(os.path.dirname(dist_file), "radius_summary")
            os.makedirs(out_dir, exist_ok=True)

            # 读取
            dist_mat = pd.read_csv(dist_file, index_col=0)
            idmap = pd.read_csv(idmap_file)

            # 检查列
            if "unique_id" not in idmap.columns or "freq" not in idmap.columns:
                print(f"[WARN] {sample_name} id_map 缺少 unique_id 或 freq 列，跳过该样本")
                continue

            # 将 freq 转 float，unique_id 保为字符串
            idmap["freq"] = pd.to_numeric(idmap["freq"], errors="coerce").fillna(0.0).astype(float)

            # 安全对齐：有些矩阵内的 id 可能在 idmap 中找不到 -> 使用 reindex，缺失 freq 设为 0 并发出警告
            seq_ids = list(dist_mat.index.astype(str))
            idmap_indexed = idmap.set_index("unique_id")
            missing = [x for x in seq_ids if x not in idmap_indexed.index]
            if missing:
                print(f"[WARN] {sample_name} 有 {len(missing)} 个矩阵 unique_id 在 id_map 中找不到；这些序列的 freq 将视为 0")
            # reindex 保留顺序
            idmap_for_mat = idmap_indexed.reindex(seq_ids)
            # 对缺失的行用 0 填充 freq
            if "freq" not in idmap_for_mat.columns:
                idmap_for_mat["freq"] = 0.0
            idmap_for_mat["freq"] = idmap_for_mat["freq"].fillna(0.0).astype(float)
            freq_map = idmap_for_mat["freq"].to_dict()

            # 计算每个 radius
            summary_records = []
            for R in radius_list:
                print(f"  radius = {R}")
                try:
                    clusters = cluster_by_radius(dist_mat, R)
                    summary = summarize_clusters(clusters, freq_map)
                    summary["radius"] = R
                    summary_records.append(summary)
                except Exception as e_inner:
                    print(f"[ERROR] sample {sample_name} radius={R} 计算失败: {e_inner}")
                    traceback.print_exc()
                    # 记录一个空/标记行，方便结果表格对齐
                    summary_records.append({
                        "n_sequence": None,
                        "n_clusters": None,
                        "avg_cluster_size": None,
                        "largest_cluster_freq": None,
                        "largest_cluster_size": None,
                        "effective_cluster_number": None,
                        "gini_index": None,
                        "top3_cluster_freq_sum": None,
                        "top10_cluster_freq_sum": None,
                        "radius": R
                    })

            # 输出 summary CSV
            out_csv = os.path.join(out_dir, f"{sample_name}_radius_summary.csv")
            pd.DataFrame(summary_records).to_csv(out_csv, index=False)
            print(f"[OK] 已生成 {out_csv}")

        except Exception as e:
            # 捕获单个样本的任意异常但继续下一个样本
            print(f"[ERROR] 处理样本 {dist_file} 时发生未捕获异常: {e}")
            traceback.print_exc()
            continue

############################################
#  主程序入口
############################################
def main():
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python TCRdist_radius_batch.py <batch_dir> <radius_list>")
        print("  radius_list 例如: 1,3,5,7,10")
        sys.exit(1)

    batch_dir = sys.argv[1]
    # parse radius list, 忽略空项，支持浮点或整数
    radius_list = []
    for x in sys.argv[2].split(","):
        xs = x.strip()
        if xs == "":
            continue
        try:
            # 尝试 float 转换再看是否是整数
            r = float(xs)
            if r.is_integer():
                r = int(r)
            radius_list.append(r)
        except:
            print(f"[WARN] 无法解析 radius 值: {xs}，将跳过")
    if not radius_list:
        print("[ERROR] 未解析到有效的 radius 列表")
        sys.exit(1)

    print(f"[INFO] 批次目录: {batch_dir}")
    print(f"[INFO] radius 列表: {radius_list}")

    process_batch(batch_dir, radius_list)
    print("\n[DONE] 所有样本处理完成")

if __name__ == "__main__":
    main()
