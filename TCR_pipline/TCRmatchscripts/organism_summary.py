# -*- coding: utf-8 -*-
import os
import pandas as pd
import argparse


## 调用 python TCR_Project/scripts/TCRmatchscripts/organism_summary.py -d /haplox/users/xuliu/TCR_Project/Results/250529PCRBias_TCR/TCRmatchResult
# ========= argparse =========
parser = argparse.ArgumentParser(description="汇总 tumor_output 中的 matched 文件")
parser.add_argument(
    "-d", "--dir",
    required=True,
    help="TCRmatchResult 目录，例如: /haplox/users/xuliu/TCR_Project/Results/250519_TCRcell/TCRmatchResult"
)
parser.add_argument(
    "-o", "--output",
    default="tumor_file_summary.csv",
    help="输出文件名 (默认: %(default)s)，会保存在 --dir 下"
)
args = parser.parse_args()

base_dir = args.dir
final_output = os.path.join(base_dir, args.output)
# ========= argparse 结束 =========

summary_records = []

# 遍历样本目录
for sample_folder in sorted(os.listdir(base_dir)):
    sample_path = os.path.join(base_dir, sample_folder)
    tumor_path = os.path.join(sample_path, "tumor_output")
    if not os.path.isdir(tumor_path):
        continue

    matched_files = [f for f in os.listdir(tumor_path) if f.endswith('.tsv') and '_matched_' in f]
    for fname in matched_files:
        file_path = os.path.join(tumor_path, fname)
        try:
            df = pd.read_csv(file_path, sep='\t')
            if not set(['sequence', 'freq']).issubset(df.columns):
                continue

            df['freq'] = pd.to_numeric(df['freq'], errors='coerce').fillna(0)
            total_freq = df['freq'].sum()
            n_seq = df['sequence'].nunique()

            # 提取 matched_ 后的部分（不含扩展名）
            target_part = fname.split('_matched_')[-1].replace('.tsv', '')

            summary_records.append({
                'sample': sample_folder,
                'target_name': target_part,
                'total_freq': total_freq,
                'n_sequence': n_seq
            })

        except Exception as e:
            print("⚠️ 文件处理失败: {} 错误: {}".format(file_path, e))

# 输出为 CSV
if summary_records:
    df_summary = pd.DataFrame(summary_records)
    df_summary = df_summary.sort_values(by=['sample', 'target_name'])
    df_summary.to_csv(final_output, index=False)
    print("✅ 汇总完成: {}".format(final_output))
else:
    print("❌ 未找到任何有效数据")
