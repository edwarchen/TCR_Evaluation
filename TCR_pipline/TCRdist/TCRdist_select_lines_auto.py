#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

ROOT_DIR = "/haplox/users/xuliu/TCR_Project/"
## 将同一批次不同文件下的TCRdist_radius_auto.py计算得到的文件选定列根据半径进行外部合并

############################################
# 提取样本名（从 <sample>_radius_summary.csv）
############################################
def extract_sample_name_from_filename(fn):
    """
    输入文件名： <sample>_radius_summary.csv
    返回 sample 名。
    """
    if fn.endswith("_radius_summary.csv"):
        return fn.replace("_radius_summary.csv", "")
    return os.path.splitext(fn)[0]


############################################
# 防止列名冲突
############################################
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


############################################
# 主程序
############################################
if len(sys.argv) < 3:
    print("Usage: python TCRdist_select_lines_auto.py <batch_dir> <column_name>")
    print("Example: python TCRdist_select_lines_auto.py /haplox/users/xuliu/TCR_Project/Results/250712_TCRcell/ top10_cluster_freq_sum")
    sys.exit(1)

batch_name = sys.argv[1]
target_col = sys.argv[2]

batch_dir = os.path.join(ROOT_DIR, batch_name)

if not os.path.isdir(batch_dir):
    print(f"Error: {batch_dir} is not a valid directory")
    sys.exit(1)

tcrdist_dir = os.path.join(batch_dir, "TCRdist")
os.makedirs(tcrdist_dir, exist_ok=True)

merged_df = None
included_samples = []
skipped_files = []
exist_names = set()


############################################
# 遍历批次目录：只匹配新的 <sample>_radius_summary.csv
############################################
for root, dirs, files in os.walk(batch_dir):
    # ⭐ 跳过旧目录 TCRdist_output/radius_summary/
    if "radius_summary" in root.split("/"):
        continue

    for fn in files:

        # ⭐ 只处理 <sample>_radius_summary.csv
        if not fn.endswith("_radius_summary.csv"):
            continue

        full = os.path.join(root, fn)

        try:
            df = pd.read_csv(full)
        except Exception as e:
            skipped_files.append((full, f"read error {e}"))
            continue

        if "radius" not in df.columns or target_col not in df.columns:
            skipped_files.append((full, f"missing col '{target_col}'"))
            continue

        df["radius"] = pd.to_numeric(df["radius"], errors="coerce")

        # ⭐ 样本名直接从文件名提取
        raw_name = extract_sample_name_from_filename(fn)
        sample_name = make_unique_name(raw_name, exist_names)

        tmp = df[["radius", target_col]].copy()
        tmp = tmp.rename(columns={target_col: sample_name})
        tmp = tmp.set_index("radius")

        if merged_df is None:
            merged_df = tmp
        else:
            merged_df = merged_df.join(tmp, how="outer")

        included_samples.append(sample_name)


############################################
# 输出结果
############################################
if merged_df is None:
    print("No valid *_radius_summary.csv found.")
    sys.exit(0)

merged_df = merged_df.sort_index().reset_index().rename(columns={"index": "radius"})

out_file = f"merged_{target_col}_by_radius.csv"
out_path = os.path.join(tcrdist_dir, out_file)
merged_df.to_csv(out_path, index=False)

print(f"\n合并完成：{out_path}")
print(f"纳入样本数量：{len(included_samples)}")

if skipped_files:
    print("\n跳过的文件：")
    for path, reason in skipped_files:
        print(f" - {path} : {reason}")
