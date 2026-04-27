#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
import os
import numpy as np
import math
from collections import defaultdict


def get_freq_stat(lst):
    total, length = np.sum(lst), len(lst)
    return total, length, total * math.log(length)


INPUT_COLS = ['freq', 'antigen.species']

vdjmatch_metadata_file = os.path.abspath(sys.argv[1])
input_dir = os.path.abspath(sys.argv[2])
output_file = os.path.abspath(sys.argv[3])

vdjmatch_metadata_tbl = pd.read_csv(vdjmatch_metadata_file, sep='\t')

print(f'Input length: {len(vdjmatch_metadata_tbl)}')

# ❗ 不再过滤 TRA/TRB（你现在是合法的）
vdjmatch_metadata_tbl.drop_duplicates(subset=['sample_id'], inplace=True)

print(f'After drop duplicates: {len(vdjmatch_metadata_tbl)}')

total_vdjmath_result_lst = []

for _, row in vdjmatch_metadata_tbl.iterrows():
    vdjmatch_result_file = os.path.join(input_dir, row['file_name'])

    if not os.path.exists(vdjmatch_result_file):
        print(f'{row["file_name"]} do not exist')
        continue

    sample_id = row['sample_id']

    # ✅ 新增：解析 chain 和真实 sample
    if sample_id.endswith('_TRA'):
        chain = 'TRA'
        sample = sample_id.replace('_TRA', '')
    elif sample_id.endswith('_TRB'):
        chain = 'TRB'
        sample = sample_id.replace('_TRB', '')
    else:
        chain = 'UNK'
        sample = sample_id

    vdjmatch_result_tbl = pd.read_csv(
        vdjmatch_result_file,
        sep='\t',
        usecols=INPUT_COLS
    )

    vdjmatch_grouped_tble = vdjmatch_result_tbl.groupby(['antigen.species']).agg({'freq': list}).reset_index()

    vdjmatch_grouped_tble['total_freq'], vdjmatch_grouped_tble['count'], vdjmatch_grouped_tble['score'] = zip(
        *vdjmatch_grouped_tble['freq'].map(get_freq_stat)
    )

    # ✅ 关键新增列
    vdjmatch_grouped_tble['sample'] = sample
    vdjmatch_grouped_tble['chain'] = chain

    vdjmatch_grouped_tble.drop(columns=['freq'], inplace=True)

    total_vdjmath_result_lst.append(vdjmatch_grouped_tble)


# 合并
total_vdjmath_result_tbl = pd.concat(total_vdjmath_result_lst)

# ✅ 输出
total_vdjmath_result_tbl.to_csv(output_file, index=None)

print(">>> Done!")