import pandas as pd, os

# 单行命令式转换
input_file = '/haplox/users/xuliu/TCR_Project/SampleList/file2022_complete.tsv'
output_dir = '/haplox/users/xuliu/TCR_Project/Results/2022_TCRcell'
output_file = f'{output_dir}/metadata.tsv'

# 确保目录存在
os.makedirs(output_dir, exist_ok=True)

# 读取、交换列、重命名、保存
df = pd.read_csv(input_file, sep='\t')
df = df.iloc[:, ::-1]  # 交换列位置
df.columns = ['sample_id', 'file_name']
df.to_csv(output_file, sep='\t', index=False)

print(f"转换完成: {output_file}")