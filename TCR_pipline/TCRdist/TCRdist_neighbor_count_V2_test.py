#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

# 指定目录
DATA_DIR = "/haplox/users/xuliu/TCR_Project/downloadfile/TCRdist_output"
OUTFILE = os.path.join(DATA_DIR, "distance_distribution_per_sample.csv")

# -------------------------------
# 查找矩阵文件
# -------------------------------
def find_matrix_files(data_dir):
    files = []
    for f in os.listdir(data_dir):
        if f.endswith("_tcrdist_matrix.csv"):
            files.append(os.path.join(data_dir, f))
    return files

# -------------------------------
# 处理单个文件
# -------------------------------
def process_one_file(fpath):
    sample = os.path.basename(fpath).replace("_tcrdist_matrix.csv", "")
    df = pd.read_csv(fpath, index_col=0)

    # 转长格式
    df_long = df.stack().reset_index()
    df_long.columns = ["seq_i", "seq_j", "dist"]

    # 去掉自己对自己的距离
    df_long = df_long[df_long["seq_i"] != df_long["seq_j"]]

    # 统计每个距离的数量
    stat = df_long["dist"].value_counts().reset_index()
    stat.columns = ["dist", "count"]
    stat["Sample"] = sample

    return stat

# -------------------------------
# 主程序
# -------------------------------
if __name__ == "__main__":
    print(f"[INFO] 扫描目录 {DATA_DIR} 下所有 *_tcrdist_matrix.csv 文件")
    files = find_matrix_files(DATA_DIR)
    print(f"[INFO] 发现 {len(files)} 个矩阵文件")

    if not files:
        print("没有发现可处理的文件，脚本结束。")
        sys.exit(0)

    all_stats = []
    for f in files:
        print(f"[INFO] 处理 {os.path.basename(f)} ...")
        res = process_one_file(f)
        all_stats.append(res)

    final_df = pd.concat(all_stats, ignore_index=True)
    final_df = final_df.sort_values(["Sample", "dist"])
    final_df.to_csv(OUTFILE, index=False)
    print(f"[OK] 已输出统计结果：{OUTFILE}")
    print(final_df.head())
