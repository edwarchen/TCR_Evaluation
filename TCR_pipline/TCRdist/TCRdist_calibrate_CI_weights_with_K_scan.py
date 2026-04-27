#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import random
import numpy as np
import pandas as pd
from collections import deque
from itertools import combinations

#python TCR_Project/scripts/TCRdist/TCRdist_calibrate_CI_weights_with_K_scan.py  --sample_list /haplox/users/xuliu/TCR_Project/SampleList/Derived_analysis.tsv --radius 30
#[INFO] Selected 15 samples: ['S084_SZ20250529017NOT-e_ntdna_genome_235619', 'S168_SZ20250409054WHB-b_gdna_genome_235222', 'S495_SZ20230818130-1_celldna_genome_136880', 
# 'S461_SZ20230818099-1_celldna_genome_136846', 'S488_SZ20230818125-1_celldna_genome_136873', 'S077_SZ20250529024PCT-d_atdna_genome_235612', 
# 'S131_SZ20250313013WHB-f_gdna_genome_235185', 'S471_SZ20230818104-8_celldna_genome_136856', 'S113_SZ20250709031NOT-1_ntdna_genome_235648', 
# 'S106_SZ20250703111TUT-f_ttdna_genome_235641', 'S486_SZ20230818124-9_celldna_genome_136871', 'S490_SZ20230818123-c_celldna_genome_136875', 
# 'S087_SZ20250529033PCT-c_atdna_genome_235622', 'S191_SZ20250416011WHB-1_gdna_genome_235245', 'S054_SZ20250123103TUT-7_ttdna_genome_230368']

EPS = 1e-8

# -------------------------
# 基于半径的簇划分
# -------------------------
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

# -------------------------
# 计算 TTT
# -------------------------
def compute_TTT(cluster_df, freq_df):
    df = cluster_df.merge(freq_df, on="unique_id", how="left")
    valid = df.groupby("cluster_id").filter(lambda x: len(x) >= 2)
    if valid.empty:
        return np.nan
    return valid["freq_norm"].sum()

# -------------------------
# 计算 DDD
# -------------------------
def compute_DDD(cluster_df, freq_df):
    df = cluster_df.merge(freq_df, on="unique_id", how="left")
    valid = df.groupby("cluster_id").filter(lambda x: len(x) >= 2)
    if valid.empty:
        return np.nan
    return valid.groupby("cluster_id")["freq_norm"].sum().max()

# -------------------------
# 计算 AAA
# -------------------------
def compute_AAA(cluster_df, dist_mat):
    valid = cluster_df.groupby("cluster_id").filter(lambda x: len(x) >= 2)
    clusters = valid["cluster_id"].unique()
    n = len(clusters)
    if n <= 1:
        return np.nan

    cluster_seqs = {c: valid.loc[valid["cluster_id"] == c, "unique_id"].values
                    for c in clusters}

    inter_dists = []
    for c1, c2 in combinations(clusters, 2):
        d = dist_mat.loc[cluster_seqs[c1], cluster_seqs[c2]].values.mean()
        inter_dists.append(d)

    inter_dist = np.mean(inter_dists)
    max_dist = dist_mat.values.max()
    return 1.0 - inter_dist / max_dist

# -------------------------
# 查找样本文件夹
# -------------------------
def find_sample_folder(root_dir, sample_id):
    for subfolder in os.listdir(root_dir):
        path = os.path.join(root_dir, subfolder, sample_id)
        if os.path.isdir(path):
            return path
    return None

# -------------------------
# 主流程
# -------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_list", required=True)
    parser.add_argument("--root_dir", default="/haplox/users/xuliu/TCR_Project/Results")
    parser.add_argument("--output_dir", default="/haplox/users/xuliu/TCR_Project/Results/TCR_CI_weight")
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--n_sample", type=int, default=15)
    parser.add_argument("--seed", type=int, default=1334)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 读取样本列表
    df_samples = pd.read_csv(args.sample_list, sep="\t")
    sample_ids = df_samples["Sample_ID"].dropna().unique().tolist()

    # 随机抽样
    random.seed(args.seed)
    selected_samples = random.sample(sample_ids, args.n_sample)
    print(f"[INFO] Selected {args.n_sample} samples: {selected_samples}")

    results_all = []

    for sid in selected_samples:
        sample_folder = find_sample_folder(args.root_dir, sid)
        if sample_folder is None:
            print(f"[WARN] Sample folder not found for {sid}, skipping.")
            continue

        dist_file = os.path.join(sample_folder, "Map_Clone_Analysis", "TCRdist_output",
                                 f"{sid}.clonotypes.TRB_pairwise_tcrdist.csv")
        freq_file = os.path.join(sample_folder, "Map_Clone_Analysis", "TCRdist_output",
                                 f"{sid}.clonotypes.TRB_unique_id_map.csv")

        if not os.path.exists(dist_file) or not os.path.exists(freq_file):
            print(f"[WARN] Required files missing for sample {sid}, skipping.")
            continue

        # 读取矩阵
        dist_mat = pd.read_csv(dist_file, index_col=0)

        # 读取频率文件并归一化
        freq_df = pd.read_csv(freq_file)
        if "freq" not in freq_df.columns or "unique_id" not in freq_df.columns:
            print(f"[ERROR] Frequency file format incorrect for sample {sid}")
            continue
        freq_df["freq_norm"] = freq_df["freq"] / freq_df["freq"].sum()

        # 聚类
        cluster_df = cluster_by_radius(dist_mat, args.radius)

        # 计算指标
        ttt = compute_TTT(cluster_df, freq_df)
        ddd = compute_DDD(cluster_df, freq_df)
        aaa = compute_AAA(cluster_df, dist_mat)

        results_all.append({
            "Sample_ID": sid,
            "radius": args.radius,
            "TTT": ttt,
            "DDD": ddd,
            "AAA": aaa,
            "n_cluster": cluster_df["cluster_id"].nunique()
        })

    # 保存每个样本指标
    result_df = pd.DataFrame(results_all)
    sample_metrics_file = os.path.join(args.output_dir, f"AAA_TTT_DDD_random{args.n_sample}_radius{args.radius}.csv")
    result_df.to_csv(sample_metrics_file, index=False)
    print(f"[DONE] Saved per-sample metrics: {sample_metrics_file}")

    # =========================
    # 计算方差和权重
    # =========================
    var_dict = result_df[["TTT", "DDD", "AAA"]].var(ddof=1)
    weights = var_dict / var_dict.sum()

    weight_df = pd.DataFrame({
        "Metric": var_dict.index,
        "Variance": var_dict.values,
        "Weight": weights.values
    })
    weight_file = os.path.join(args.output_dir, f"TCR_CI_weights_radius{args.radius}.csv")
    weight_df.to_csv(weight_file, index=False)
    print(f"[DONE] Saved final weights: {weight_file}")
    print(weight_df)


if __name__ == "__main__":
    main()
