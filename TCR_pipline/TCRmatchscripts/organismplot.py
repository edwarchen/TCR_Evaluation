# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 防止无图形界面报错
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# 调用 python TCR_Project/scripts/TCRmatchscripts/organismplot.py -d /haplox/users/xuliu/TCR_Project/Results/250529PCRBias_TCR/TCRmatchResult/

# ========= 新增 argparse =========
parser = argparse.ArgumentParser(description="统计 matched TSV 并绘图")
parser.add_argument(
    "-d", "--dir",
    required=True,
    help="TCRmatchResult 目录，例如: /haplox/users/xuliu/TCR_Project/Results/250519_TCRcell/TCRmatchResult"
)
parser.add_argument(
    "-o", "--output",
    default="tumor_output",
    help="匹配结果子文件夹 (默认: %(default)s)"
)
parser.add_argument(
    "-p", "--plots",
    default="plots",
    help="绘图结果子文件夹 (默认: %(default)s)"
)
args = parser.parse_args()

base_dir = args.dir
output_subfolder = args.output
plot_subfolder = args.plots
# ========= argparse 结束 =========

# 遍历每个样本目录
for sample_folder in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample_folder)
    output_path = os.path.join(sample_path, output_subfolder)
    if not os.path.isdir(output_path):
        continue

    matched_files = [f for f in os.listdir(output_path) if f.endswith('.tsv') and '_matched_' in f]
    if not matched_files:
        print("跳过无tsv文件: {}".format(output_path))
        continue

    records = []
    for f in matched_files:
        file_path = os.path.join(output_path, f)
        try:
            df = pd.read_csv(file_path, sep='\t')
            if not {'sequence', 'freq'}.issubset(df.columns):
                continue

            target = f.split("_matched_")[-1].replace(".tsv", "")
            df['freq'] = pd.to_numeric(df['freq'], errors='coerce').fillna(0)

            total_freq = df['freq'].sum()
            n_sequence = df['sequence'].nunique()

            records.append({
                "target_name": target,
                "total_freq": total_freq,
                "n_sequence": n_sequence
            })
        except Exception as e:
            print("⚠️ 文件读取失败: {} 错误: {}".format(file_path, e))

    if not records:
        print("⚠️ 无可用数据: {}".format(sample_folder))
        continue

    summary_df = pd.DataFrame(records).sort_values(by='total_freq', ascending=False)

    # 创建图输出路径
    plot_dir = os.path.join(output_path, plot_subfolder)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # 保存表格
    summary_csv = os.path.join(plot_dir, "target_summary_freq_and_sequence.csv")
    summary_df.to_csv(summary_csv, index=False)

    # === 图1: total_freq 柱状图 ===
    plt.figure(figsize=(10, max(4, len(summary_df) * 0.4)))
    sns.barplot(data=summary_df, x="total_freq", y="target_name", palette="Blues_d")
    plt.xlabel("Total Freq")
    plt.ylabel("Target Name")
    plt.title("Sample: {} | Total Freq".format(sample_folder))
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "{}_total_freq.png".format(sample_folder)), bbox_inches='tight')
    plt.close()

    # === 图2: n_sequence 柱状图 ===
    summary_df_seq = summary_df.sort_values(by='n_sequence', ascending=False)
    plt.figure(figsize=(10, max(4, len(summary_df_seq) * 0.4)))
    sns.barplot(data=summary_df_seq, x="n_sequence", y="target_name", palette="Greens_d")
    plt.xlabel("Number of Unique Sequences")
    plt.ylabel("Target Name")
    plt.title("Sample: {} | Unique Sequence Count".format(sample_folder))
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "{}_n_sequence.png".format(sample_folder)), bbox_inches='tight')
    plt.close()

    print("✅ {} 绘图完成: 共 {} 项".format(sample_folder, len(summary_df)))
