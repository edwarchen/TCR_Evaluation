# -*- coding: utf-8 -*-
import pandas as pd

# 定义文件路径
files = ["/haplox/users/xuliu/TCR_Project/Results/250519_TCRcell/TCRmatchResult/final_combined.tsv", "/haplox/users/xuliu/TCR_Project/Results/250303_bdszB_TCR/TCRmatchResult/final_combined.tsv"]

# 用于存储处理后的DataFrame
filtered_dfs = []

for file in files:
    # 读取文件
    df = pd.read_csv(file, sep='\t')
    
    # 获取第一列名
    first_col = df.columns[0]
    
    # 找出包含关键字的列
    keyword_cols = [col for col in df.columns if "Mycobacterium tuberculosis" in col or "Lung" in col]
    
    # 构造最终列列表（保持顺序）
    selected_cols = [first_col] + keyword_cols

    # 选择这些列
    filtered_df = df[selected_cols]
    filtered_dfs.append(filtered_df)

# 合并为一个DataFrame
combined_df = pd.concat(filtered_dfs, ignore_index=True)

# 保存为新的TSV文件
combined_df.to_csv("/haplox/users/xuliu/TCR_Project/Results/Lung/combined_filtered.tsv", sep='\t', index=False)

print("处理完成：已生成 combined_filtered.tsv")
