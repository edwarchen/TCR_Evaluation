# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')  # 设置非交互式后端
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# 读取数据
file_path = "/haplox/users/xuliu/TCR_Project/Results/250519_TCRcell/TCRmatchResult/final_combined.tsv"
df = pd.read_csv(file_path, sep='\t')

# 数据预处理
def melt_data(df, cols, suffix):
    melted = pd.melt(df, id_vars=["sample"], value_vars=cols, var_name="Category", value_name=suffix)
    melted["Category"] = melted["Category"].str.replace(r"_{}$".format(suffix), "", regex=True)
    return melted

count_cols = [col for col in df.columns if col.endswith("_count")]
freq_cols = [col for col in df.columns if col.endswith("_freq")]
count_df = melt_data(df, count_cols, "count")
freq_df = melt_data(df, freq_cols, "freq")
merged_df = pd.merge(count_df, freq_df, on=["sample", "Category"])

# 绘图函数
def plot_data(data, y_col, title):
    plt.figure(figsize=(80, 12))
    
    samples = data["sample"].unique()
    colors = sns.color_palette("hls", n_colors=len(samples))
    sample_palette = dict(zip(samples, colors))

    sns.barplot(
        data=data,
        x="Category",
        y=y_col,
        hue="sample",
        palette=sample_palette,
        dodge=True
    )

    
    plt.title("{} by Category and Sample".format(y_col.capitalize()), fontsize=14)
    plt.xlabel("Category", fontsize=12)
    plt.ylabel(y_col.capitalize(), fontsize=12)
    plt.xticks(rotation=90, fontsize=8)
    plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.gca().ticklabel_format(axis="y", style="plain")

    plt.tight_layout()
    plt.savefig("/haplox/users/xuliu/TCR_Project/Results/250519_TCRcell/TCRmatchResult/{}.png".format(title), dpi=150, bbox_inches="tight")
    plt.close()


# 生成图表
plot_data(merged_df, "count", "count_plot")
plot_data(merged_df, "freq", "freq_plot")