#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd

MATRIX_ROOT = "/haplox/users/xuliu/TCR_Project/Results"
target_seq = "CASSSYNEQFF"
RADIUS = [2, 4, 6, 8, 10]
OUTPUT = "tcr_neighbor_summary.csv"


def find_matrix_files(root):
    files = []
    for dp, dn, fn in os.walk(root):
        for f in fn:
            if f.endswith("_tcrdist_matrix.csv") or f.endswith("_cdr3only_levenshtein.csv"):
                files.append(os.path.join(dp, f))
    return files


def analyze_file(fpath, target_seq, radii):
    try:
        df = pd.read_csv(fpath, index_col=0)
    except:
        print(f"[WARN] 无法读取：{fpath}")
        return None

    seqs = df.index.tolist()

    if target_seq not in seqs:
        return None

    distances = df.loc[target_seq, :]

    # -----------【关键修复】-----------
    if not isinstance(distances, pd.Series):
        print(f"[WARN] {fpath} → 发现重复的 CDR3 行名或格式异常，跳过该样本")
        return None
    # -----------------------------------

    distances = pd.to_numeric(distances, errors="coerce").fillna(999999)

    counts = {}
    for r in radii:
        counts[r] = (distances <= r).sum() - 1  # 去掉自己
    return counts


if __name__ == "__main__":
    files = find_matrix_files(MATRIX_ROOT)
    print(f"[INFO] 找到 {len(files)} 个矩阵文件")

    results = []

    for f in files:
        sample = os.path.basename(f).replace("_tcrdist_matrix.csv", "").replace("_cdr3only_levenshtein.csv","")
        res = analyze_file(f, target_seq, RADIUS)
        if res:
            row = {"Sample": sample}
            row.update(res)
            results.append(row)

    out_df = pd.DataFrame(results)
    out_df.to_csv(OUTPUT, index=False)

    print(f"[OK] 已生成结果文件: {OUTPUT}")
    print(out_df)
