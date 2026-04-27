#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

ROOT_DIR = "/haplox/users/xuliu/TCR_Project/"

if len(sys.argv) < 3:
    print("Usage: python TCRdist_select_lines.py <batch_name> <column_name>")
    print("Example: python TCRdist_select_lines.py 250712_TCRcell top10_cluster_freq_sum")
    sys.exit(1)

batch_name = sys.argv[1]
target_col = sys.argv[2]   # ⭐ 你要提取的列名

batch_dir = os.path.join(ROOT_DIR, batch_name)

if not os.path.isdir(batch_dir):
    print(f"Error: {batch_dir} is not a valid directory")
    sys.exit(1)

merged_df = None
included_samples = []
skipped_files = []

# 防止列名冲突
def make_unique_name(name, exist_set):
    if name not in exist_set:
        exist_set.add(name)
        return name
    i = 1
    while f"{name}_{i}" in exist_set:
        i += 1
    new = f"{name}_{i}"
    exist_set.add(new)
    return new

exist_names = set()

for root, dirs, files in os.walk(batch_dir):
    for fn in files:
        if not fn.endswith("_radius_summary.csv"):
            continue
        full = os.path.join(root, fn)

        try:
            df = pd.read_csv(full)
        except Exception as e:
            skipped_files.append((full, f"read error {e}"))
            continue

        # 必须要有 radius 和 你指定的目标列
        if "radius" not in df.columns or target_col not in df.columns:
            skipped_files.append((full, f"missing column '{target_col}'"))
            continue

        try:
            df["radius"] = pd.to_numeric(df["radius"], errors="coerce")
        except:
            pass

        base = fn.replace("_radius_summary.csv", "")
        sample_name = make_unique_name(base, exist_names)

        tmp = df[["radius", target_col]].copy()
        tmp = tmp.rename(columns={target_col: sample_name})
        tmp = tmp.set_index("radius")

        if merged_df is None:
            merged_df = tmp
        else:
            merged_df = merged_df.join(tmp, how="outer")

        included_samples.append(sample_name)

if merged_df is None:
    print("[WARN] No valid summary files found.")
    sys.exit(0)

merged_df = (
    merged_df
    .sort_index()
    .reset_index()
    .rename(columns={"index": "radius"})
)

out_file = f"merged_{target_col}_by_radius.csv"
out_path = os.path.join(batch_dir, out_file)
merged_df.to_csv(out_path, index=False)

print(f"\n合并完成: {out_path}")
print(f"纳入样本数量: {len(included_samples)}")

if skipped_files:
    print("跳过的文件:")
    for path, reason in skipped_files:
        print(f" - {path} : {reason}")
