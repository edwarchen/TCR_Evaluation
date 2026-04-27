#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
统计所有样本的 TCR 序列（cdr3_b_aa）出现频次，找出 public TCR。
扫描所有 *_tcrdist_matrix.csv 和 *_cdr3only_levenshtein.csv 文件。
输出：/haplox/users/xuliu/TCR_Project/Results/public_TCR_top.csv
"""

import os
import pandas as pd
from collections import defaultdict

# -----------------------------------------------------------
# 目录设置
# -----------------------------------------------------------
ROOT_DIR = "/haplox/users/xuliu/TCR_Project/Results"
OUTPUT_FILE = os.path.join(ROOT_DIR, "public_TCR_top.csv")

# -----------------------------------------------------------
# 扫描矩阵文件
# -----------------------------------------------------------
matrix_files = []
for dirpath, dirnames, filenames in os.walk(ROOT_DIR):
    for fn in filenames:
        if fn.endswith("_tcrdist_matrix.csv") or fn.endswith("_cdr3only_levenshtein.csv"):
            matrix_files.append(os.path.join(dirpath, fn))

print(f"[INFO] 共找到 {len(matrix_files)} 个矩阵文件")

# -----------------------------------------------------------
# 统计序列出现频率
# -----------------------------------------------------------
tcr_to_samples = defaultdict(set)

for f in matrix_files:
    sample_name = os.path.basename(f).split("_cdr3")[0].split("_tcrdist")[0]

    try:
        df = pd.read_csv(f, index_col=0)
    except Exception as e:
        print(f"[WARN] 无法读取文件 {f}, 跳过。错误: {e}")
        continue

    # 行名 = CDR3 序列
    seqs = df.index.tolist()

    for s in seqs:
        tcr_to_samples[s].add(sample_name)

print(f"[INFO] 共汇总到 {len(tcr_to_samples)} 条独特序列")

# -----------------------------------------------------------
# 输出结果整理
# -----------------------------------------------------------
rows = []
for seq, sset in tcr_to_samples.items():
    rows.append({
        "CDR3": seq,
        "n_samples": len(sset),
        "samples": "|".join(sorted(sset))
    })

out_df = pd.DataFrame(rows)
out_df = out_df.sort_values("n_samples", ascending=False)

out_df.to_csv(OUTPUT_FILE, index=False)
print(f"[OK] 已保存 public TCR 统计结果至: {OUTPUT_FILE}")
