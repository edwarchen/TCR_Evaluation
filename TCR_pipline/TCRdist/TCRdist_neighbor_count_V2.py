#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

ROOT = "/haplox/users/xuliu/TCR_Project/Results"
OUTFILE = "distance_distribution_per_sample.csv"

#python distance_dist_stat.py
#python distance_dist_stat.py Batch_A Batch_B

def find_matrix_files(root, batch_list=None):
    """
    batch_list: 如果为 None，则扫描整个 ROOT；
                 如果是 ['20240601_batch1', '20240601_batch2'] 则只扫描对应子目录。
    """
    out = []

    # 没指定批次 → 扫描全部目录
    if not batch_list:
        for dp, dn, fn in os.walk(root):
            for f in fn:
                if f.endswith("_cdr3only_levenshtein.csv"):
                    out.append(os.path.join(dp, f))
        return out

    # 指定批次 → 只扫描对应目录
    for batch in batch_list:
        batch_path = os.path.join(root, batch)
        if not os.path.exists(batch_path):
            print(f"[WARN] 批次不存在: {batch}")
            continue

        for dp, dn, fn in os.walk(batch_path):
            for f in fn:
                if f.endswith("_cdr3only_levenshtein.csv"):
                    out.append(os.path.join(dp, f))

    return out


def process_one_file(fpath):
    sample = os.path.basename(fpath).replace("_cdr3only_levenshtein.csv", "")

    df = pd.read_csv(fpath, index_col=0)

    df_long = df.stack().reset_index()
    df_long.columns = ["seq_i", "seq_j", "dist"]

    df_long = df_long[df_long["seq_i"] != df_long["seq_j"]]

    stat = df_long["dist"].value_counts().reset_index()
    stat.columns = ["dist", "count"]
    stat["Sample"] = sample

    return stat


if __name__ == "__main__":

    # -------------------------------
    # 参数解析
    # -------------------------------
    # 用法示例：
    #   python script.py batch1 batch2
    #
    # 无参数 → 扫描整个 ROOT
    # -------------------------------
    batch_list = sys.argv[1:] if len(sys.argv) > 1 else None
    if batch_list:
        print(f"[INFO] 仅统计批次: {batch_list}")
    else:
        print("[INFO] 未指定批次，将统计整个 ROOT")

    files = find_matrix_files(ROOT, batch_list=batch_list)
    print(f"[INFO] 发现 {len(files)} 个矩阵文件")

    if not files:
        print("没有发现可处理的文件，脚本结束。")
        sys.exit(0)

    all_stats = []

    for f in files:
        print(f"[INFO] 处理 {f} ...")
        res = process_one_file(f)
        all_stats.append(res)

    final_df = pd.concat(all_stats, ignore_index=True)
    final_df = final_df.sort_values(["Sample", "dist"])

    final_df.to_csv(OUTFILE, index=False)
    print(f"[OK] 已输出统计结果：{OUTFILE}")
    print(final_df.head())
