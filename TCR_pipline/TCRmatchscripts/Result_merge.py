#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import pandas as pd
import os
import sys

# 工作目录作为第一个命令行参数
if len(sys.argv) < 2:
    print("Usage: python merge_results.py <working_dir>")
    sys.exit(1)

wd = sys.argv[1]

# 输入文件路径
combined_stat_path = os.path.join(wd, "all_combined_stat.csv")
final_combined_path = os.path.join(wd, "TCRmatchResult/final_combined.tsv")

# 输出文件路径
output_path = os.path.join(wd, "Results_vdj_tcrmatch.tsv")

# 读取 CSV / TSV 文件
combined_stat = pd.read_csv(combined_stat_path, encoding='utf-8')
final_combined = pd.read_csv(final_combined_path, sep="\t", encoding='utf-8')

# 检查列名
if "Sample_ID" not in combined_stat.columns:
    raise ValueError("all_combined_stat.csv 缺少列 Sample_ID")
if "sample" not in final_combined.columns:
    raise ValueError("final_combined.tsv 缺少列 sample")

# 合并
merged = pd.merge(combined_stat, final_combined, left_on="Sample_ID", right_on="sample", how="inner")

# 保存结果
merged.to_csv(output_path, sep="\t", index=False)
print("Merged result saved to: {}".format(output_path))