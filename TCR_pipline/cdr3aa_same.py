#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import pandas as pd
import itertools

root_dir = "/haplox/users/xuliu/TCR_Project/TCR_BCR_RD_Proj/TR_6sample_2"
topN_list = [100, 500, 1000]
freq_cutoff = 0

# --------- 函数 ---------
def parse_sample_info(sample_name):
    return {
        "method": "Self" if "Self_Devel" in sample_name else "Takara",
        "sample": "Samp1",
        "rep": "Rep1" if "Rep1" in sample_name else "Rep2"
    }

def load_top_df(file, topN):
    df = pd.read_csv(file, sep="\t")

    df = df.dropna(subset=["cdr3aa"]).drop_duplicates("cdr3aa")

    if "freq" in df.columns:
        df = df.sort_values("freq", ascending=False)
    elif "count" in df.columns:
        df = df.sort_values("count", ascending=False)
        df["freq"] = df["count"] / df["count"].sum()
    else:
        raise ValueError(f"{file} 没有freq或count列")

    df = df[df["freq"] > freq_cutoff]

    return df.head(topN)[["cdr3aa", "freq"]]

# --------- 收集文件 ---------
sample_files = {}

for sample in os.listdir(root_dir):
    map_dir = os.path.join(root_dir, sample, "Map_Clone_Analysis")
    if not os.path.isdir(map_dir):
        continue

    files = glob.glob(os.path.join(map_dir, "convert.*.txt"))

    for f in files:
        if "TRA" in f:
            chain = "TRA"
        elif "TRB" in f:
            chain = "TRB"
        else:
            continue

        sample_files.setdefault(chain, {})[sample] = f

# --------- 主计算 ---------
overlap_results = []
all_freq_all = []

for chain in ["TRA", "TRB"]:

    if chain not in sample_files:
        continue

    samples = list(sample_files[chain].keys())

    # ⭐ 自动两两组合
    for g1, g2 in itertools.combinations(samples, 2):

        file1 = sample_files[chain][g1]
        file2 = sample_files[chain][g2]

        info1 = parse_sample_info(g1)
        info2 = parse_sample_info(g2)

        for topN in topN_list:

            df1 = load_top_df(file1, topN)
            df2 = load_top_df(file2, topN)

            set1 = set(df1["cdr3aa"])
            set2 = set(df2["cdr3aa"])

            shared = set1 & set2

            # --------- overlap指标 ---------
            overlap_count = len(shared)

            # ⭐ 更合理（推荐）
            jaccard = overlap_count / len(set1 | set2) if len(set1 | set2) > 0 else 0

            overlap_results.append({
                "chain": chain,
                "group1": g1,
                "group2": g2,
                "method1": info1["method"],
                "method2": info2["method"],
                "rep1": info1["rep"],
                "rep2": info2["rep"],
                "topN": topN,
                "overlap_count": overlap_count,
                "jaccard_index": jaccard
            })

            # --------- freq ---------
            df1 = df1.rename(columns={"freq": "freq_g1"})
            df2 = df2.rename(columns={"freq": "freq_g2"})

            merged = pd.merge(df1, df2, on="cdr3aa", how="outer")

            merged["freq_g1"] = merged["freq_g1"].fillna(0)
            merged["freq_g2"] = merged["freq_g2"].fillna(0)

            def classify(row):
                if row["freq_g1"] > 0 and row["freq_g2"] > 0:
                    return "shared"
                elif row["freq_g1"] > 0:
                    return "only_g1"
                else:
                    return "only_g2"

            merged["type"] = merged.apply(classify, axis=1)

            merged["chain"] = chain
            merged["group1"] = g1
            merged["group2"] = g2
            merged["topN"] = topN

            all_freq_all.append(merged)

# --------- 输出 ---------
df_overlap = pd.DataFrame(overlap_results)
df_overlap.to_csv("cdr3_overlap_count.tsv", sep="\t", index=False)

if len(all_freq_all) > 0:
    df_all_freq = pd.concat(all_freq_all)
    df_all_freq.to_csv("cdr3_all_freq.tsv", sep="\t", index=False)

print("✅ 完成（升级版）")