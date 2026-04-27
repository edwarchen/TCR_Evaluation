#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import networkx as nx
import os
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description="Convert wide TCRdist matrix to Cytoscape edge, node & cluster files"
    )
    parser.add_argument(
        "--matrix",
        required=True,
        help="pairwise TCRdist matrix (rows & columns are node IDs)"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="output directory"
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=12,
        help="TCRdist threshold (default: 12)"
    )
    parser.add_argument(
        "--sep",
        default=",",
        help="file separator (default: ,)"
    )

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1. 读取宽格式矩阵
    mat = pd.read_csv(args.matrix, sep=args.sep, index_col=0)

    # 基本健壮性检查
    assert (mat.columns == mat.index).all(), \
        "Row names and column names must be identical"

    nodes = mat.index.tolist()

    edges = []

    # 2. 只遍历上三角，避免重复
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            d = mat.iat[i, j]
            if pd.isna(d):
                continue
            if d <= args.radius:
                edges.append({
                    "source": nodes[i],
                    "target": nodes[j],
                    "weight": float(d)
                })

    edge_df = pd.DataFrame(edges)

    edge_file = os.path.join(
        args.outdir, f"edges_R{args.radius}.tsv"
    )
    edge_df.to_csv(edge_file, sep="\t", index=False)

    # 3. 用边构图，算连通分量（singleton 会形成 size=1 的 component）
    G = nx.Graph()
    G.add_nodes_from(nodes)

    for _, row in edge_df.iterrows():
        G.add_edge(row["source"], row["target"])

    clusters = []
    singleton_nodes = set()

    for i, comp in enumerate(nx.connected_components(G), start=1):
        comp = list(comp)
        if len(comp) == 1:
            singleton_nodes.add(comp[0])
        for node in comp:
            clusters.append({
                "node_id": node,
                "cluster_id": f"cluster_{i}"
            })

    cluster_df = pd.DataFrame(clusters)

    cluster_file = os.path.join(
        args.outdir, f"clusters_R{args.radius}.tsv"
    )
    cluster_df.to_csv(cluster_file, sep="\t", index=False)

    # ===== NEW: nodes.tsv（两列，Cytoscape 可直接导入）=====
    nodes_df = pd.DataFrame({
        "node_id": nodes,
        "is_singleton": [node in singleton_nodes for node in nodes]
    })

    node_file = os.path.join(args.outdir, "nodes.tsv")
    nodes_df.to_csv(node_file, sep="\t", index=False)

    print("Done.")
    print(f"Nodes   : {node_file}")
    print(f"Edges   : {edge_file}")
    print(f"Clusters: {cluster_file}")
    print(f"Total nodes     : {len(nodes)}")
    print(f"Total edges     : {len(edge_df)}")
    print(f"Total clusters  : {cluster_df['cluster_id'].nunique()}")
    print(f"Singleton nodes : {len(singleton_nodes)}")

if __name__ == "__main__":
    main()
