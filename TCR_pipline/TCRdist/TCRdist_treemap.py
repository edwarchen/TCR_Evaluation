#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import squarify
import matplotlib.colors as mcolors


# -----------------------------
# 参数解析
# -----------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dist_matrix", required=True)
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out_prefix", required=True)
    return parser.parse_args()


# -----------------------------
# 读取距离矩阵
# -----------------------------
def read_distance_matrix(path):
    if path.endswith(".csv"):
        df = pd.read_csv(path, index_col=0)
    else:
        df = pd.read_csv(path, sep="\t", index_col=0)

    df = df.loc[df.index, df.index]
    return df


# -----------------------------
# 半径聚类
# -----------------------------
def radius_clustering(dist_df, radius):
    G = nx.Graph()
    ids = dist_df.index.tolist()
    G.add_nodes_from(ids)

    mat = dist_df.values
    n = len(ids)

    for i in range(n):
        for j in range(i + 1, n):
            if mat[i, j] <= radius:
                G.add_edge(ids[i], ids[j])

    clusters = list(nx.connected_components(G))

    cluster_map = {}
    for idx, comp in enumerate(clusters):
        cid = f"C{idx + 1}"
        for uid in comp:
            cluster_map[uid] = cid

    return pd.Series(cluster_map, name="cluster_id")


# -----------------------------
# 主流程
# -----------------------------
def main():
    args = parse_args()

    print("▶ Reading distance matrix")
    dist_df = read_distance_matrix(args.dist_matrix)

    print("▶ Reading metadata")
    meta = pd.read_csv(args.metadata)
    meta = meta[["unique_id", "freq"]]

    print(f"▶ Radius clustering (r = {args.radius})")
    cluster_series = radius_clustering(dist_df, args.radius)

    df = (
        cluster_series
        .reset_index()
        .rename(columns={"index": "unique_id"})
        .merge(meta, on="unique_id", how="left")
    )

    cluster_df = (
        df.groupby("cluster_id")
        .agg(
            n_seq=("unique_id", "size"),
            cluster_freq=("freq", "sum")
        )
        .reset_index()
        .sort_values("cluster_freq", ascending=False)
    )

    # 仅取前500
    cluster_df = cluster_df.head(500).copy()

    total_freq = cluster_df["cluster_freq"].sum()
    cluster_df["relative_freq"] = cluster_df["cluster_freq"] / total_freq

    # singleton 不显示数字
    cluster_df["label"] = cluster_df["n_seq"].apply(
        lambda x: str(x) if x > 1 else ""
    )

    cluster_df.to_csv(
        f"{args.out_prefix}/top500_cluster_summary.tsv",
        sep="\t",
        index=False
    )

    print("▶ Top 500 cluster summary saved")

    # -----------------------------
    # 🎨 手动马卡龙配色 + 透明度递减
    # -----------------------------
    base_colors = [
        "#F28E8E",  # 柔红
        "#F5B971",  # 柔橙
        "#F9E07F",  # 奶黄
        "#8ED1A6",  # 薄荷绿
        "#77DCDD",  # 青蓝
        "#6B98C4",  # 深蓝
        "#C99BFF",  # 紫
        "#E07BB5",  # 洋红
        "#BABABA",  # 中性灰
        "white",  # 柔棕
    ]

    colors = []
    num_clusters = len(cluster_df)

    alpha_start = 0.95
    alpha_decay = 0.15
    min_alpha = 0.35

    for i in range(num_clusters):
        cycle = i // 10
        position = i % 10

        alpha = max(alpha_start - cycle * alpha_decay, min_alpha)

        rgb = mcolors.to_rgb(base_colors[position])
        rgba = (*rgb, alpha)

        colors.append(rgba)

    # -----------------------------
    # Treemap
    # -----------------------------
    plt.figure(figsize=(12, 10))

    squarify.plot(
        sizes=cluster_df["relative_freq"],
        label=cluster_df["label"],
        color=colors,
        edgecolor="white",
        linewidth=1.0,
        text_kwargs={
            "fontsize": 10,
            "color": "black"
        }
    )

    plt.axis("off")
    plt.title(
        f"TCRdist radius = {args.radius}\nTop 500 clusters (relative frequency)",
        fontsize=14
    )

    plt.tight_layout()
    plt.savefig(
        f"{args.out_prefix}/top500_treemap.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    print("▶ Treemap saved (manual pastel + alpha decay)")


if __name__ == "__main__":
    main()
