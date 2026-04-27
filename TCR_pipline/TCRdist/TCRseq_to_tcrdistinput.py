#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import re

"""
批次运行示例：
python TCRseq_to_tcrdistinput.py 250716_TCRcell_supple

功能说明：
1. 自动进入 /haplox/users/xuliu/TCR_Project/Results/<批次号>/
2. 自动发现所有 convert.*.clonotypes.TRB.txt 文件
3. 读取 clonetype 文件（列名：count freq cdr3nt cdr3aa v d j VEnd DStart DEnd JStart）
4. 自动将 v/j 补全等位基因，如 TRBV12-3 → TRBV12-3*01
5. 输出同目录下 *_formatted.tsv
"""

ROOT_DIR = "/haplox/users/xuliu/TCR_Project/Results/"


# -----------------------------
# 修正规因子格式 TRBV12-3 → TRBV12-3*01
# -----------------------------
def fix_allele(gene):
    if pd.isna(gene):
        return gene

    gene = gene.strip()

    if "*" in gene:
        return gene   # 已经是等位基因

    m = re.match(r"(TR[AB][VDJ][0-9\-]+)", gene)
    if m:
        return m.group(1) + "*01"

    return gene


# -----------------------------
# 处理单个 clonotype 文件
# -----------------------------
def process_clonotype_file(infile):
    df = pd.read_csv(infile, sep="\t", dtype=str)

    # 你的 clonetype 文件真实列名如下，必须完全匹配
    required = ["count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j",
                "VEnd", "DStart", "DEnd", "JStart"]

    missing = [x for x in required if x not in df.columns]
    if missing:
        print(f"[WARN] {os.path.basename(infile)} 缺少列：{missing}，跳过。")
        return None

    # 补写 V/J 基因等位基因
    df["v"] = df["v"].apply(fix_allele)
    df["j"] = df["j"].apply(fix_allele)

    # 输出给 tcrdist 的字段
    out = df[["cdr3aa", "cdr3nt", "v", "j", "count", "freq"]].copy()

    # 转数值
    out["count"] = pd.to_numeric(out["count"], errors="coerce").fillna(0)
    out["freq"] = pd.to_numeric(out["freq"], errors="coerce").fillna(0)

    return out


# -----------------------------
# 批量操作
# -----------------------------
def process_batch(batch_name):
    batch_path = os.path.join(ROOT_DIR, batch_name)

    if not os.path.isdir(batch_path):
        print(f"[ERROR] 批次目录不存在：{batch_path}")
        return

    print(f"[INFO] 扫描批次目录：{batch_path}")

    for root, dirs, files in os.walk(batch_path):
        for f in files:
            if f.startswith("convert.") and f.endswith(".clonotypes.TRB.txt"):
                infile = os.path.join(root, f)
                outfile = infile.replace(".txt", "_formatted.tsv")

                print(f"[INFO] 处理文件：{infile}")

                df2 = process_clonotype_file(infile)
                if df2 is None:
                    continue

                df2.to_csv(outfile, sep="\t", index=False)
                print(f"[OK] 输出文件：{outfile}")


# -----------------------------
# 主程序
# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("使用方法：python TCRseq_to_tcrdistinput.py <批次号>")
        sys.exit(1)

    batch_name = sys.argv[1]
    process_batch(batch_name)
