#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import pandas as pd

#将samplelist.txt对应的样本名中的cdr3aa筛选count≥10，C开头F结尾，不包含_或*的cdr3aa行并重新计算筛选后的freq
# 配置路径
results_dir = "/x03_haplox/users/xuliu/TCR_Project/Results"
sample_list_file = "/x03_haplox/users/xuliu/TCR_Project/olga_prob/samplelist.txt"
output_folder_name = "Filter_cdr3aa"  # 所有输出统一保存到这里

# 输出文件夹路径
output_dir = os.path.join(results_dir, output_folder_name)

# 如果存在就先删除
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)

# 创建新的输出文件夹
os.makedirs(output_dir, exist_ok=True)

# 读取样本列表
with open(sample_list_file) as f:
    samples_to_process = [line.strip().split()[0] for line in f.readlines()]

# 遍历 Results 下所有文件夹
for root_folder in os.listdir(results_dir):
    folder_path = os.path.join(results_dir, root_folder)
    if not os.path.isdir(folder_path):
        continue

    # 遍历文件夹下的子文件夹，查找样本名称
    for sample_folder in os.listdir(folder_path):
        if sample_folder not in samples_to_process:
            continue
        sample_path = os.path.join(folder_path, sample_folder)
        map_clone_path = os.path.join(sample_path, "Map_Clone_Analysis")
        if not os.path.exists(map_clone_path):
            continue

        # 找到 convert.*.clonotypes.TRB.txt 文件
        for f in os.listdir(map_clone_path):
            if f.startswith("convert.") and f.endswith(".clonotypes.TRB.txt"):
                file_path = os.path.join(map_clone_path, f)
                df = pd.read_csv(file_path, sep="\t")

                # 筛选条件
                df_filtered = df[
                    (df['count'] >= 10) &
                    df['cdr3aa'].str.startswith('C') &
                    df['cdr3aa'].str.endswith('F') &
                    (~df['cdr3aa'].str.contains(r'[_*]'))
                ].copy()

                # 重新计算 freq，存入新列 cdr3aa_freq
                total_freq = df_filtered['freq'].sum()
                if total_freq > 0:
                    df_filtered['cdr3aa_freq'] = df_filtered['freq'] / total_freq
                else:
                    df_filtered['cdr3aa_freq'] = 0

                # 输出文件到统一文件夹
                output_file = os.path.join(output_dir, f"filtered.{sample_folder}.{f}")
                df_filtered.to_csv(output_file, sep="\t", index=False)
                print(f"Processed {file_path} -> {output_file}")
