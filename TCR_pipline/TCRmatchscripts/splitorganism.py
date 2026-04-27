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
    if os.path.isdir(sample_path):
        input_file = os.path.join(sample_path, 'merge_renamed.tsv')
        if os.path.exists(input_file):
            try:
                df = pd.read_csv(input_file, sep='\t')
                required_cols = ['trimmed_input_sequence', 'freq', 'renamed_organism']
                if not all(col in df.columns for col in required_cols):
                    print("缺少列: {}".format(input_file))
                    continue
                if df.empty:
                    print("空文件: {}".format(input_file))
                    continue

                # 拆分 renamed_organism
                expanded_rows = []
                for row in df.itertuples(index=False):
                    organism_field = getattr(row, 'renamed_organism')
                    if pd.isnull(organism_field) or organism_field.strip() == '':
                        continue
                    organisms = str(organism_field).split(',')
                    for org in organisms:
                        expanded_rows.append({
                            'input': getattr(row, 'trimmed_input_sequence'),
                            'freq': getattr(row, 'freq'),
                            'organism': org.strip()
                        })

                df_expanded = pd.DataFrame(expanded_rows)

                # 去重：每个 input-organism 保留一次
                df_unique = df_expanded.drop_duplicates(subset=['input', 'organism'])

                # 添加 count 列：统计每个 input-freq-organism 出现的次数（每条唯一序列算一次）
                df_counted = (
                    df_unique.groupby(['input', 'freq', 'organism'])
                    .size()
                    .reset_index(name='count')
                )

                # 排序
                df_counted = df_counted.sort_values(by='freq', ascending=False)

                # 输出
                output_file = os.path.join(sample_path, 'processed.tsv')
                df_counted.to_csv(output_file, sep='\t', index=False)
                print("完成样本: {}".format(sample_folder))

            except Exception as e:
                print("处理失败: {} 错误: {}".format(input_file, e))
