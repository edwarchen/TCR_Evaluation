# -*- coding: utf-8 -*-
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--base_dir', required=True, help='Base directory for TCR results')
args = parser.parse_args()

base_dir = args.base_dir

for sample_folder in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample_folder)
    processed_file = os.path.join(sample_path, 'processed.tsv')
    if os.path.exists(processed_file):
        try:
            df = pd.read_csv(processed_file, sep='\t')
            if df.empty:
                print("空文件: {}".format(processed_file))
                continue

            # 删除 input 列
            df = df.drop(columns=['input'])

            # 转换 freq 和 count 为数字类型
            df['freq'] = pd.to_numeric(df['freq'], errors='coerce').fillna(0)
            df['count'] = pd.to_numeric(df['count'], errors='coerce').fillna(0)

            # 按 organism 汇总
            summary = (
                df.groupby('organism', as_index=False)
                .agg({'freq': 'sum', 'count': 'sum'})
            )

            # 输出 summary 文件
            output_file = os.path.join(sample_path, 'organism_summary.tsv')
            summary.to_csv(output_file, sep='\t', index=False)
            print("完成汇总: {}".format(sample_folder))

        except Exception as e:
            print("统计失败: {} 错误: {}".format(processed_file, e))
