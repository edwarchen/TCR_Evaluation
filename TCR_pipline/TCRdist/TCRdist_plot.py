#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

# ======== 参数区 ========
input_csv = "/haplox/users/xuliu/TCR_Project/Results/250716_TCRcell/S058_SZ20250529010TUT-6_ttdna_genome_235593/Map_Clone_Analysis/TCRdist_output/S058_SZ20250529010TUT-6_ttdna_genome_235593_cdr3only_levenshtein_5000.csv"
sample_id = "S058_SZ20250529010TUT-6"
output_prefix = input_csv.replace(".csv", "")

# ======== Step 1: 读取距离矩阵 ========
print(f"[INFO] 读取 {input_csv}")
dist_df = pd.read_csv(input_csv, index_col=0)
print(f"[INFO] 距离矩阵形状: {dist_df.shape}")

# ======== Step 2: UMAP降维 ========
print("[INFO] 运行 UMAP 降维...")
embedding = umap.UMAP(
    n_neighbors=15,     # 可调，越大越关注全局结构
    min_dist=0.5,       # 越小越成团
    metric='precomputed',
    random_state=42
).fit_transform(dist_df.values)

umap_df = pd.DataFrame(embedding, columns=["UMAP1", "UMAP2"], index=dist_df.index)

# ======== Step 3: 基于DBSCAN聚类 ========
print("[INFO] 基于 DBSCAN 进行聚类...")
clustering = DBSCAN(eps=2.5, min_samples=5, metric='precomputed').fit(dist_df.values)
umap_df["Cluster"] = clustering.labels_

# ======== Step 4: 绘制可视化 ========
plt.figure(figsize=(8,6))
unique_clusters = np.unique(umap_df["Cluster"])
colors = plt.cm.get_cmap("tab20", len(unique_clusters))

for i, cluster_id in enumerate(unique_clusters):
    subset = umap_df[umap_df["Cluster"] == cluster_id]
    label = f"Cluster {cluster_id}" if cluster_id != -1 else "Noise"
    plt.scatter(subset["UMAP1"], subset["UMAP2"], s=10, alpha=0.7, label=label, color=colors(i))

plt.title(f"{sample_id} TCR Similarity Clustering (UMAP + DBSCAN)")
plt.legend(markerscale=2, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.tight_layout()
plt.savefig(f"{output_prefix}_umap_dbscan_dist05_15.png", dpi=300)
plt.show()

# ======== Step 5: 导出结果 ========
out_csv = f"{output_prefix}_umap_cluster.csv"
umap_df.to_csv(out_csv)
print(f"[OK] 输出聚类坐标表: {out_csv}")
