#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: python combine_all_stat2.py <input_dir> <output_file>")
    sys.exit(1)

input_dir = os.path.abspath(sys.argv[1])
output_file = os.path.abspath(sys.argv[2])

DEFAULT_ID = 'Sample_ID'

file_dict = {
    'fastp_qc.csv': {'sample':'Sample', 'chain':None},
    'vdjtools_stat/diversity.strict.resampled.txt': {'sample':'sample_id', 'chain':'chain'},
    'total_result.tsv': {'sample':'sample', 'chain':'chain'},
    'cdr3aa_stat_summary.csv': {'sample':'Sample_ID', 'chain':None},
    'enriched_cdr3aa.csv': {'sample':'sample_id', 'chain':None},
    'aa_stat.csv': {'sample':'Sample_ID', 'chain':'Chain'},
    'length_stat.csv': {'sample':'Sample_ID', 'chain':'Chain'},
    'convergence_stat.csv': {'sample':'Sample_ID', 'chain':'Chain'},
    'total_V.csv': {'sample':'Sample_ID', 'chain':'chain'},
    'total_J.csv': {'sample':'Sample_ID', 'chain':'chain'},
    'reshaped_vdjdb_anno_summary.csv': {'sample':'sample', 'chain':'chain'},
}

sample_chain_map = {}
total_table_lst = []

# ==========================
# 处理主文件
# ==========================
for file_name, cols in file_dict.items():

    if file_name == 'fastp_qc.csv':
        continue

    file_path = os.path.join(input_dir, file_name)
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} does not exist, skipping.")
        continue

    sep = '\t' if file_name.endswith(('.tsv','.txt')) else ','
    df = pd.read_csv(file_path, sep=sep)

    # ========= 自动识别 sample 列 =========
    sample_col = cols['sample']
    if sample_col not in df.columns:
        for alt in ['Sample_ID', 'sample_id', 'Sample', 'sample']:
            if alt in df.columns:
                sample_col = alt
                break
        else:
            raise KeyError(f"No sample column found in {file_name}")

    # ========= 统一 Sample_ID =========
    if sample_col != DEFAULT_ID:
        df[DEFAULT_ID] = df[sample_col]
        df.drop(columns=[sample_col], inplace=True)

    # ========= 特殊拆分 =========
    if file_name in ['cdr3aa_stat_summary.csv', 'enriched_cdr3aa.csv']:

        orig_id = df[DEFAULT_ID].copy()

        tmp = df[DEFAULT_ID].astype(str).str.extract(r'^(.*)_(TRA|TRB)$')

        df[DEFAULT_ID] = tmp[0]
        df['chain'] = tmp[1]

        df[DEFAULT_ID] = df[DEFAULT_ID].fillna(orig_id)
        df['chain'] = df['chain'].fillna('unknown')

    else:
        df[DEFAULT_ID] = df[DEFAULT_ID].astype(str).str.replace(r'_(TRA|TRB)$', '', regex=True)

    # ========= 统一 chain =========
    if cols['chain'] is not None and cols['chain'] in df.columns:
        df['chain'] = df[cols['chain']]
        if cols['chain'] != 'chain':
            df.drop(columns=[cols['chain']], inplace=True)

    if 'chain' not in df.columns:
        df['chain'] = 'unknown'

    # ========= 记录 sample-chain =========
    for s, c in zip(df[DEFAULT_ID], df['chain']):
        sample_chain_map.setdefault(s, set()).add(c)

    total_table_lst.append(df)

# ==========================
# 处理 fastp
# ==========================
fastp_file = os.path.join(input_dir, 'fastp_qc.csv')
fastp_table = pd.DataFrame()

if os.path.exists(fastp_file):
    df_fastp = pd.read_csv(fastp_file)

    df_fastp[DEFAULT_ID] = df_fastp['Sample']
    df_fastp.drop(columns=['Sample'], inplace=True)

    df_fastp[DEFAULT_ID] = df_fastp[DEFAULT_ID].str.replace(r'_(TRA|TRB)$', '', regex=True)

    fastp_rows = []
    for _, row in df_fastp.iterrows():
        s = row[DEFAULT_ID]
        chains = sample_chain_map.get(s, ['unknown'])
        for ch in chains:
            new_row = row.copy()
            new_row['chain'] = ch
            fastp_rows.append(new_row)

    fastp_table = pd.DataFrame(fastp_rows)
else:
    print("Warning: fastp_qc.csv not found.")

# ==========================
# 合并
# ==========================
all_tables = [fastp_table] + total_table_lst if not fastp_table.empty else total_table_lst
total_table = pd.concat(all_tables, ignore_index=True)

total_table = total_table.groupby([DEFAULT_ID, 'chain'], as_index=False).first()

# ==========================
# fastq比例
# ==========================
if 'Clean_Reads_Num(M)' in total_table.columns and 'total_clone_reads' in total_table.columns:
    total_table['Clean_Reads'] = total_table['Clean_Reads_Num(M)'] * 1e6
    total_table['chain_ratio_to_fastq'] = total_table.apply(
        lambda x: x['total_clone_reads'] / x['Clean_Reads'] if x['Clean_Reads'] > 0 else 0,
        axis=1
    )

# ==========================
# TRA/TRB比例
# ==========================
if 'total_clone_reads' in total_table.columns:

    ratio_df = total_table[total_table['chain'].isin(['TRA', 'TRB'])].copy()

    ratio_df['total_reads_per_sample'] = ratio_df.groupby('Sample_ID')['total_clone_reads'].transform('sum')

    ratio_df['chain_ratio_in_sample'] = ratio_df.apply(
        lambda x: x['total_clone_reads'] / x['total_reads_per_sample']
        if x['total_reads_per_sample'] > 0 else 0,
        axis=1
    )

    ratio_df = ratio_df[['Sample_ID', 'chain', 'chain_ratio_in_sample']]

    total_table = total_table.merge(
        ratio_df,
        on=['Sample_ID', 'chain'],
        how='left'
    )

else:
    print("Warning: total_clone_reads not found.")

# ==========================
# 排序列
# ==========================
cols = total_table.columns.tolist()

front = [c for c in ['Sample_ID', 'chain'] if c in cols]
others = [c for c in cols if c not in front]

total_table = total_table[front + others]

# ==========================
# 输出
# ==========================
os.makedirs(os.path.dirname(output_file), exist_ok=True)
total_table.to_csv(output_file, sep=",", index=False)

print("✅ Done:", output_file)