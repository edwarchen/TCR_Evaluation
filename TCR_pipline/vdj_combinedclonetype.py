#!/usr/bin/env python3
# coding: utf-8

#调用：python vdj_combinedclonetype.py E20250718（批次号）  得到：同一路径下的convert*.cdr3aa_freq.tsv


import os
import sys
import pandas as pd

# ✅ 写死的根路径
PROJECT_ROOT = "/haplox/users/xuliu/TCR_Project/Results"

def find_target_files(batch_path):
    """
    查找 batch 目录下每个 sample 的 Map_Clone_Analysis 目录中的 convert*.clonotypes.TRB.txt 文件
    """
    target_files = []
    for sample in os.listdir(batch_path):
        sample_path = os.path.join(batch_path, sample)
        if not os.path.isdir(sample_path):
            continue
        map_clone_dir = os.path.join(sample_path, "Map_Clone_Analysis")
        if not os.path.isdir(map_clone_dir):
            continue
        for fname in os.listdir(map_clone_dir):
            if fname.startswith("convert") and fname.endswith(".clonotypes.TRB.txt"):
                target_files.append(os.path.join(map_clone_dir, fname))
    return target_files

def process_file(file_path):
    """
    读取并聚合 cdr3aa-freq，输出 .cdr3aa_freq.tsv
    """
    print(f"Processing: {file_path}")
    try:
        df = pd.read_csv(file_path, sep='\t', dtype={'cdr3aa': str})
        if 'cdr3aa' in df.columns and 'freq' in df.columns:
            df_agg = df.groupby('cdr3aa', as_index=False)['freq'].sum()
            out_file = file_path.replace(".clonotypes.TRB.txt", ".cdr3aa_freq.tsv")
            df_agg.to_csv(out_file, sep='\t', index=False)
            print(f"Saved to: {out_file}")
        else:
            print(f"Skipped: missing 'cdr3aa' or 'freq' in {file_path}")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python aggregate_cdr3aa_freq.py <batch_name>")
        sys.exit(1)

    batch_name = sys.argv[1]
    batch_path = os.path.join(PROJECT_ROOT, batch_name)

    if not os.path.isdir(batch_path):
        print(f"Error: batch directory not found: {batch_path}")
        sys.exit(1)

    files = find_target_files(batch_path)
    if not files:
        print("No target files found.")
        return

    for f in files:
        process_file(f)

if __name__ == "__main__":
    main()
