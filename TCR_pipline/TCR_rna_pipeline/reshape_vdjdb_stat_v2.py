#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys
import os
import numpy as np
import math

# ===== 统计函数 =====
def get_freq_stat(lst):
    total, length = np.sum(lst), len(lst)
    if length == 0:
        return 0, 0, 0
    return total, length, total * math.log(length)


# ===== 参数 =====
SELECTED_ANTIGEN_SPECIES = [
    "HIV","YFV","CMV","EBV","HCV",
    "HomoSapiens","InfluenzaA","SARS-CoV-2","TriticumAestivum"
]

input_file = os.path.abspath(sys.argv[1])
output_file = os.path.abspath(sys.argv[2])

# ===== 读取数据 =====
input_table = pd.read_csv(input_file)

# ===== 检查 chain 列 =====
if "chain" not in input_table.columns:
    raise ValueError("❌ 输入文件中没有 chain 列，请检查上一步输出")

# ===== 过滤 antigen =====
input_table = input_table[
    input_table['antigen.species'].isin(SELECTED_ANTIGEN_SPECIES)
]

# ===== 🚀 关键修改：加入 chain =====
reshaped_table = input_table.pivot_table(
    index=['sample', 'chain'],
    columns='antigen.species',
    values=['total_freq', 'count'],
    fill_value=0
)

# ===== 展平成普通列 =====
reshaped_table = reshaped_table.reset_index()

# ===== 扁平化列名 =====
reshaped_table.columns = [
    '_'.join(col).strip('_') if isinstance(col, tuple) else col
    for col in reshaped_table.columns
]

# ===== 保存 =====
reshaped_table.to_csv(output_file, index=False)

print(">>> Done!")