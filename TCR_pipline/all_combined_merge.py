# -*- coding: utf-8 -*-
import pandas as pd
import os

files = [
    "/haplox/users/xuliu/TCR_Project/Results/250519_TCRcell/all_combined_stat.csv",
    "/haplox/users/xuliu/TCR_Project/Results/250303_bdszB_TCR/all_combined_stat.csv"
]

# 读取第一个文件，保留列顺序
df0 = pd.read_csv(files[0])
columns_order = df0.columns.tolist()

# 读取第二个文件，强制列顺序匹配第一个
df1 = pd.read_csv(files[1])
df1 = df1[columns_order]  # 重新排序列

# 合并
combined_df = pd.concat([df0, df1], ignore_index=True)

# 保存
combined_df.to_csv("/haplox/users/xuliu/TCR_Project/Results/20250605/all_combined_stat_merged.csv", index=False, encoding='utf-8')

print("✅ 合并完成，列顺序保持一致。")