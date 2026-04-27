# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

base_dir = "/haplox/users/xuliu/TCR_Project/Results/250519_TCRcell/TCRmatchResult"
output_subfolder = "tumor_output"
plot_subfolder = "plots"
plot_filename = "combined_targets_freq.png"
plot_field = "freq"  # 只绘制 freq

for sample_folder in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample_folder)
    output_path = os.path.join(sample_path, output_subfolder)
    if not os.path.isdir(output_path):
        continue

    matched_files = [f for f in os.listdir(output_path) if f.endswith('.tsv') and '_matched_' in f]
    if not matched_files:
        print("跳过无符合条件的tsv文件夹: {}".format(output_path))
        continue

    combined_df = pd.DataFrame()

    for f in matched_files:
        file_path = os.path.join(output_path, f)
        try:
            df = pd.read_csv(file_path, sep='\t')
            if 'sequence' in df.columns and plot_field in df.columns:
                target = f.split("_matched_")[-1].replace(".tsv", "")
                df['target_rename'] = target
                combined_df = pd.concat([combined_df, df], ignore_index=True)
        except Exception as e:
            print("⚠️ 文件读取失败: {} 错误: {}".format(file_path, e))

    if combined_df.empty:
        print("无数据可绘图: {}".format(sample_folder))
        continue

    # 生成唯一序列名
    def make_unique_seqs(df):
        seen = {}
        new_seqs = []
        for seq in df['sequence']:
            if seq not in seen:
                seen[seq] = 0
                new_seqs.append(seq)
            else:
                seen[seq] += 1
                new_seqs.append(seq + '+' * seen[seq])
        return new_seqs

    combined_df['sequence_unique'] = None
    for target, group in combined_df.groupby('target_rename'):
        combined_df.loc[group.index, 'sequence_unique'] = make_unique_seqs(group)

    # 去掉 + 用于标红判断
    combined_df['sequence_clean'] = combined_df['sequence_unique'].str.replace(r'\++$', '', regex=True)

    # 分组数据用于绘图
    grouped_df = combined_df[['target_rename', 'sequence_unique', 'sequence_clean', 'freq']].copy()
    grouped_df = grouped_df.sort_values(['target_rename', plot_field], ascending=[True, False])

    grouped_df['seq_label'] = grouped_df['target_rename'] + ' | ' + grouped_df['sequence_unique']
    grouped_df['seq_label'] = pd.Categorical(grouped_df['seq_label'], categories=grouped_df['seq_label'], ordered=True)

    # 🔴 标红逻辑：只根据 sequence_clean 是否出现在多个 target_rename 中
    red_seqs = (
        combined_df.groupby('sequence_clean')['target_rename']
        .nunique()
        .reset_index(name='n_targets')
    )
    red_keys = set(red_seqs[red_seqs['n_targets'] > 1]['sequence_clean'])

    # 输出标红序列
    plot_dir = os.path.join(output_path, plot_subfolder)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    red_out_path = os.path.join(plot_dir, "red_marked_sequences.csv")
    pd.DataFrame(sorted(red_keys), columns=["sequence_clean"]).to_csv(red_out_path, index=False)

    # Y轴标签颜色
    yticklabels = []
    for _, row in grouped_df.iterrows():
        seq_clean = str(row['sequence_clean'])
        seq_label = row['seq_label']
        color = 'red' if seq_clean in red_keys else 'black'
        yticklabels.append((seq_label, color))

    # 配色
    unique_targets = grouped_df['target_rename'].unique()
    palette = sns.color_palette("tab10", n_colors=len(unique_targets))
    color_map = dict(zip(unique_targets, palette))

    # 绘图
    plt.figure(figsize=(12, max(4, len(grouped_df)*0.25)))
    sns.set(style="whitegrid")
    sns.barplot(
        x=plot_field,
        y="seq_label",
        hue="target_rename",
        data=grouped_df,
        dodge=False,
        palette=color_map
    )

    ax = plt.gca()
    ax.set_yticks(range(len(yticklabels)))
    ax.set_yticklabels([lbl for lbl, _ in yticklabels])
    for ticklabel, color in zip(ax.get_yticklabels(), [c for _, c in yticklabels]):
        ticklabel.set_color(color)
        ticklabel.set_fontsize(8)

    plt.title("Sample: {} (Plot: {})".format(sample_folder, plot_field))
    plt.xlabel(plot_field.capitalize())
    plt.ylabel("Target | Sequence")

    plt.legend(
        title="target_organism",
        bbox_to_anchor=(1.01, 1),
        loc='upper left',
        fontsize=8,
        borderaxespad=0.
    )
    plt.tight_layout()

    out_path = os.path.join(plot_dir, plot_filename)
    plt.savefig(out_path, bbox_inches='tight')
    plt.close()

    print("✅ 绘图完成: {}".format(out_path))
