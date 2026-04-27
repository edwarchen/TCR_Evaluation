# -*- coding: utf-8 -*-
import os
import pandas as pd
import re

# 要处理的样本文件夹路径
sample_path = "/haplox/users/xuliu/TCR_Project/Results/250121_TCRcell/TCRmatchResult/S082_SZ20250110101-1_celldna_genome_219456"

# 匹配的关键字
target_organism = "Influenza A"

def safe_filename(name):
    name = name.lower()
    name = re.sub(r'\W+', '_', name)
    return name.strip('_')

output_suffix = safe_filename(target_organism) + "_hits.tsv"


input_file = os.path.join(sample_path, 'merge_renamed.tsv')

if os.path.exists(input_file):
    try:
        df = pd.read_csv(input_file, sep='\t')
        if df.empty:
            print("空文件: {}".format(input_file))
        else:
            required_cols = ['renamed_organism', 'freq', 'trimmed_input_sequence']
            missing = [col for col in required_cols if col not in df.columns]
            if missing:
                print("缺少列 {} in {}".format(', '.join(missing), input_file))
            else:
                matched_rows = []
                for i, row in df.iterrows():
                    organisms = str(row['renamed_organism']).split(',')
                    organisms = [o.strip() for o in organisms]
                    if target_organism in organisms:
                        matched_rows.append({
                            'sequence': row['trimmed_input_sequence'],
                            'freq': row['freq']
                        })

                if matched_rows:
                    df_matched = pd.DataFrame(matched_rows)
                    df_matched = df_matched.drop_duplicates(subset=['sequence', 'freq'])
                    output_file = os.path.join(sample_path, output_suffix)
                    df_matched.to_csv(output_file, sep='\t', index=False)
                    print("✅ 匹配并去重完成: {}".format(output_file))
                else:
                    print("ℹ️ 无匹配: {}".format(input_file))
    except Exception as e:
        print("❌ 处理失败 {}: {}".format(input_file, e))
else:
    print("❌ 未找到文件: {}".format(input_file))
