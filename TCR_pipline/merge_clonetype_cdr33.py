#!/usr/bin/env python3
# coding: utf-8

# 调用方式：python merge_clonetype_cdr33.py 20250718_TCRcell[batch_id] 20250718_allclonetype.tsv

import os
import sys
import pandas as pd
from glob import glob

# ✅ 写死根路径
PROJECT_ROOT = "/haplox/users/xuliu/TCR_Project/Results"

def collect_files(batch):
    result_files = []
    batch_path = os.path.join(PROJECT_ROOT, batch)
    if not os.path.isdir(batch_path):
        print(f"[Error] Batch not found: {batch_path}")
        sys.exit(1)

    sample_dirs = [os.path.join(batch_path, d) for d in os.listdir(batch_path)
                   if os.path.isdir(os.path.join(batch_path, d))]
    for sample_path in sample_dirs:
        map_clone_dir = os.path.join(sample_path, "Map_Clone_Analysis")
        if not os.path.isdir(map_clone_dir):
            continue
        for f in glob(os.path.join(map_clone_dir, "convert*.cdr3aa_freq.tsv")):
            result_files.append((batch, os.path.basename(sample_path), f))
    return result_files, batch_path

def merge_tables(file_info_list):
    all_dfs = []
    for batch, sample, file_path in file_info_list:
        try:
            df = pd.read_csv(file_path, sep='\t', dtype={'cdr3aa': str})
            if 'cdr3aa' in df.columns and 'freq' in df.columns:
                df['batch'] = batch
                df['sample'] = sample
                all_dfs.append(df[['batch', 'sample', 'cdr3aa', 'freq']])
            else:
                print(f"[Skip] Missing columns in: {file_path}")
        except Exception as e:
            print(f"[Error] Failed reading {file_path}: {e}")
    if all_dfs:
        return pd.concat(all_dfs, ignore_index=True)
    else:
        return pd.DataFrame(columns=['batch', 'sample', 'cdr3aa', 'freq'])

def main():
    if len(sys.argv) != 3:
        print("Usage: python merge_clonetype_cdr33.py <batch> <output_filename>")
        sys.exit(1)

    batch = sys.argv[1]
    output_name = sys.argv[2]

    file_info_list, batch_path = collect_files(batch)
    print(f"Found {len(file_info_list)} result files in batch {batch}.")

    merged_df = merge_tables(file_info_list)
    if merged_df.empty:
        print("No valid data found. Exiting.")
        return

    outfile = os.path.join(batch_path, output_name)
    merged_df.to_csv(outfile, sep='\t', index=False)
    print(f"[Done] Merged result saved to: {outfile}")

if __name__ == "__main__":
    main()
