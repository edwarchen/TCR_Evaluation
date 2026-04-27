#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import networkx as nx
import os

def main():
    parser = argparse.ArgumentParser(
        description="Convert wide TCRdist matrix to Cytoscape edge, node & cluster files with freq"
    )
    parser.add_argument("--matrix", required=True,
                        help="pairwise TCRdist matrix (rows & columns are node IDs)")
    parser.add_argument("--freq_file", required=True,
                        help="Node frequency file (columns: unique_id / freq)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--radius", type=float, default=12, help="TCRdist threshold (default: 12)")
    parser.add_argument("--sep", default=",", help="File separator (default: ,)")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. 读取宽格式矩阵
    mat = pd.read_csv(args.matrix, sep=args.sep, index_col=0)
    assert (mat.columns == mat.index).all(), "Row names and column names must be identical"

    nodes = mat.index.tolist()

    # 读取 freq 文件
    freq_df = pd.read_csv(args.freq_file, sep=args.sep)
    if 'unique_id' not in freq_df.columns or 'freq' not in freq_df.columns:
        raise ValueError("freq_file must have columns: unique_id,freq")

    # ===== 输出 nodes.tsv =====
    nodes_df = pd.DataFrame({"unique_id": nodes})
    nodes_df = nodes_df.merge(freq_df[['unique_id','freq']], on='unique_id', how='left')

    # 初始化 is_singleton 列为 True，后续有边的改为 False
    nodes_df['is_singleton'] = True

    node_file = os.path.join(args.outdir, "nodes.tsv")
    nodes_df.to_csv(node_file, sep="\t", index=False)

    # ===== 构建 edges =====
    edges = []
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            d = mat.iat[i, j]
            if pd.isna(d):
                continue
            if d <= args.radius:
                edges.append({
                    "source": nodes[i],
                    "target": nodes[j],
                    "weight": float(d)
                })
                # 只要有边，节点就不是 singleton
                nodes_df.loc[nodes_df['unique_id'].isin([nodes[i], nodes[j]]), 'is_singleton'] = False

    edge_df = pd.DataFrame(edges)
    edge_file = os.path.join(args.outdir, f"edges_R{args.radius}.tsv")
    edge_df.to_csv(edge_file, sep="\t", index=False)

    # ===== 用边构图，算 clusters =====
    G = nx.Graph()
    G.add_nodes_from(nodes)
    for _, row in edge_df.iterrows():
        G.add_edge(row['source'], row['target'])

    clusters = []
    for i, comp in enumerate(nx.connected_components(G), start=1):
        for node in comp:
            clusters.append({
                "unique_id": node,
                "cluster_id": f"cluster_{i}"
            })

    cluster_df = pd.DataFrame(clusters)
    cluster_file = os.path.join(args.outdir, f"clusters_R{args.radius}.tsv")
    cluster_df.to_csv(cluster_file, sep="\t", index=False)

    # ===== 更新 nodes.tsv 的 cluster_id =====
    nodes_df = nodes_df.merge(cluster_df, on='unique_id', how='left')
    nodes_df.to_csv(node_file, sep="\t", index=False)

    print("Done.")
    print(f"Nodes   : {node_file}")
    print(f"Edges   : {edge_file}")
    print(f"Clusters: {cluster_file}")
    print(f"Total nodes   : {len(nodes_df)}")
    print(f"Total edges   : {len(edge_df)}")
    print(f"Total clusters: {cluster_df['cluster_id'].nunique()}")

if __name__ == "__main__":
    main()
