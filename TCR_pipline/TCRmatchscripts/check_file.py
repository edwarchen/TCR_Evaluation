import pandas as pd

#比较在不同cutoff时，输出的tcrmatch文件内容是否一致，小的cutoff是否包含大的cutoff的结果
# 文件路径
file_a = "/haplox/users/xuliu/TCR_Project/Results/250817_TCRcell/TCRmatchResult/S059_SZ20250704102TUT-9_ttdna_genome_238349/S059_SZ20250704102TUT-9_ttdna_genome_238349_tcrmatch.tsv"  # A 文件
file_b = "/haplox/users/xuliu/TCR_Project/Results/250817_TCRcell/TCRmatchResult_97/S059_SZ20250704102TUT-9_ttdna_genome_238349/S059_SZ20250704102TUT-9_ttdna_genome_238349_tcrmatch.tsv"  # B 文件

# 读取文件
df_a = pd.read_csv(file_a, sep="\t")
df_b = pd.read_csv(file_b, sep="\t")

# 确保只用关心的列
cols = ["trimmed_input_sequence", "match_sequence", "score"]
df_a = df_a[cols]
df_b = df_b[cols]

# 1. 检查 B 中在 A 中是否存在相同的 trimmed_input_sequence + match_sequence
merged = df_b.merge(df_a, on=["trimmed_input_sequence","match_sequence"], how='left', suffixes=('_B','_A'), indicator=True)

# 2. B中在A中没有对应sequence的行
b_not_in_a_seq = merged[merged['_merge'] == 'left_only']

# 3. B中对应sequence存在，但score不同的行
score_diff = merged[(merged['_merge'] == 'both') & (merged['score_B'] != merged['score_A'])]

# 4. A中多出来的行（B中没有对应sequence）
a_not_in_b = df_a.merge(df_b, on=["trimmed_input_sequence","match_sequence"], how='left', indicator=True)
a_not_in_b = a_not_in_b[a_not_in_b['_merge'] == 'left_only']

# 输出结果
b_not_in_a_seq.to_csv("B_not_in_A_by_sequence.tsv", sep="\t", index=False)
score_diff.to_csv("B_score_diff_in_A.tsv", sep="\t", index=False)
a_not_in_b.to_csv("A_extra.tsv", sep="\t", index=False)

print("对比完成：")
print(f"B中sequence不在A中: {len(b_not_in_a_seq)} 行")
print(f"B中score不同: {len(score_diff)} 行")
print(f"A中多出来的行: {len(a_not_in_b)} 行")
