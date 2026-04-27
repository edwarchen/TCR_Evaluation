#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Network-based clustering of TCR sequences using pairwise distance.
Author: LiuX + GPT5 (Enhanced Visualization + Cluster Size Output)
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# ---------------- 参数区 ----------------
INPUT_FILE = "/haplox/users/xuliu/TCR_Project/Results/250716_TCRcell/S058_SZ20250529010TUT-6_ttdna_genome_235593/Map_Clone_Analysis/TCRdist_output/S058_SZ20250529010TUT-6_ttdna_genome_235593_cdr3only_levenshtein_5000.csv"
DIST_THRESHOLD = 6   # <= 此距离内的序列视为相似，可根据经验调整
OUTPUT_PREFIX = INPUT_FILE.replace(".csv", "_network")

# --------------------------------------

print(f"[INFO] 读取距离矩阵: {INPUT_FILE}")
dist_df = pd.read_csv(INPUT_FILE, index_col=0)
print(f"[INFO] 矩阵大小: {dist_df.shape}")

# 删除空行或无效序列
dist_df = dist_df.dropna(how='all').fillna(9999)

# 创建 NetworkX 图
G = nx.Graph()
for seq in dist_df.index:
    G.add_node(seq)

print(f"[INFO] 按阈值 {DIST_THRESHOLD} 添加边...")
for i, seq1 in enumerate(dist_df.index):
    for j, seq2 in enumerate(dist_df.columns[i+1:], i+1):
        dist = dist_df.iloc[i, j]
        if dist <= DIST_THRESHOLD:
            G.add_edge(seq1, seq2, weight=dist)

print(f"[INFO] 图节点数: {G.number_of_nodes()}, 边数: {G.number_of_edges()}")

# 计算连通分量（每个分量视为一个簇）
components = list(nx.connected_components(G))
print(f"[INFO] 检测到 {len(components)} 个簇")

cluster_map = {}
for i, comp in enumerate(components):
    for node in comp:
        cluster_map[node] = i
nx.set_node_attributes(G, cluster_map, "Cluster")

# 计算度数最高的节点作为 anchor，并记录 cluster size
anchor_nodes = []
for i, comp in enumerate(components):
    subG = G.subgraph(comp)
    if len(subG) == 0:
        continue
    anchor = max(subG.degree, key=lambda x: x[1])[0]
    cluster_size = len(subG.nodes())
    anchor_nodes.append((i, anchor, cluster_size))

# 保存包含 cluster size 的 CSV
anchor_df = pd.DataFrame(anchor_nodes, columns=["Cluster", "Anchor", "Cluster_Size"])
anchor_df.to_csv(f"{OUTPUT_PREFIX}_anchors.csv", index=False)
print(f"[OK] Anchor + Cluster size 已保存至: {OUTPUT_PREFIX}_anchors.csv")

# ============ 🎨 可视化部分 ============
print("[INFO] 正在绘制网络...")

pos = nx.spring_layout(G, seed=42, k=0.25)

unique_clusters = sorted(set(cluster_map.values()))
cluster_colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(unique_clusters)))
node_colors = [cluster_colors[cluster_map[n] % len(cluster_colors)] for n in G.nodes()]

anchor_nodes_list = [a[1] for a in anchor_nodes]
node_sizes = [800 if n in anchor_nodes_list else 120 for n in G.nodes()]

plt.figure(figsize=(12, 10))
nx.draw_networkx_edges(G, pos, alpha=0.15, width=0.3, edge_color="gray")
nx.draw_networkx_nodes(G, pos,
                       node_color=node_colors,
                       node_size=node_sizes,
                       alpha=0.9,
                       linewidths=0.4,
                       edgecolors="black")

nx.draw_networkx_nodes(G, pos,
                       nodelist=anchor_nodes_list,
                       node_color='none',
                       node_size=[s*1.6 for s in node_sizes if s == 800],
                       edgecolors='red',
                       linewidths=1.5)

plt.title(f"TCR Similarity Network (≤ {DIST_THRESHOLD}), Clusters: {len(components)}", fontsize=13)
plt.axis("off")
plt.tight_layout()
plt.savefig(f"{OUTPUT_PREFIX}_network.png", dpi=400)
plt.close()
print(f"[DONE] 网络图已保存: {OUTPUT_PREFIX}_network.png")
