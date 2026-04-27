#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 合并TCR所有链的质控信息
import os
import sys
import glob
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: python bcr_qc.py <input_dir> <output_file>")
    sys.exit(1)

input_dir = os.path.abspath(sys.argv[1])
output_file = sys.argv[2]  # 直接指定输出文件名，不拼接目录

# -----------------------------
# 1. 递归扫描所有 .result.index.txt 文件
# -----------------------------
all_files = glob.glob(os.path.join(input_dir, "**", "*.result.index.txt"), recursive=True)

# -----------------------------
# 2. 排除链文件（.TRB. / .TRA.）
# -----------------------------
result_files = [
    f for f in all_files
    if all(x not in os.path.basename(f) for x in [".TRA.", ".TRB."])
]

if not result_files:
    print("No sample result.index.txt files found!")
    sys.exit(1)

# -----------------------------
# 3. 读取并合并所有样本文件
# -----------------------------
all_dfs = []
for file in result_files:
    df = pd.read_csv(file, sep="\t")
    
    # 确保 sample 列存在
    if "sample" not in df.columns:
        df["sample"] = os.path.basename(file).replace(".result.index.txt","")
    
    all_dfs.append(df)

merged_df = pd.concat(all_dfs, ignore_index=True)

# -----------------------------
# 4. 输出整合表格
# -----------------------------
merged_df.to_csv(output_file, sep="\t", index=False)

print("Finished. Results saved in:", output_file)