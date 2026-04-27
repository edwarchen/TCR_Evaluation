# -*- coding: utf-8 -*-
import os
import pandas as pd
import re
import sys
import argparse

reload(sys)
sys.setdefaultencoding('utf-8')

#调用 python TCR_Project/scripts/TCRmatchscripts/get_organism_all.py -d /haplox/users/xuliu/TCR_Project/Results/250529PCRBias_TCR

target_renames = ['Tumor', 'Breast cancer', 'Carcinoma', 'Cervical cancer', 'Clear cell renal carcinoma', 'Colorectal cancer','Neoantigen',
                 'Gastroesophageal cancer', 'Head and neck squamous cell carcinoma', 'Leukemia', 'Lung cancer', 'Lymphoma','Tumor associated antigen',
                 'Melanoma', 'Merkel cell carcinoma', 'Ovarian cancer', 'Pancreatic cancer', 'Oncoprotein', 'Carcinoma of uterine tube',
                 'Bile duct cancer', 'Bladder cancer', 'Endometrial cancer', 'Mycobacterium tuberculosis', 'Chronic Obstructive Lung Disease','others_cancers']

## base_dir = "/haplox/users/xuliu/TCR_Project/Results/250519_TCRcell"
## tcrmatch_dir = os.path.join(base_dir, "TCRmatchResult")
## mapping_file = "/haplox/users/xuliu/TCR_Project/scripts/TCRmatchscripts/category.tsv"
## output_subfolder = "tumor_output"

parser = argparse.ArgumentParser(description="TCRmatch结果处理")
parser.add_argument(
    "-d", "--dir", 
    required=True, 
    help="结果目录，例如 /haplox/users/xuliu/TCR_Project/Results/250519_TCRcell"
)
parser.add_argument(
    "-m", "--mapping", 
    default="/haplox/users/xuliu/TCR_Project/scripts/TCRmatchscripts/category.tsv",
    help="category 映射文件 (默认: %(default)s)"
)
parser.add_argument(
    "-o", "--output", 
    default="tumor_output",
    help="输出子文件夹 (默认: %(default)s)"
)

args = parser.parse_args()

base_dir = args.dir
tcrmatch_dir = os.path.join(base_dir, "TCRmatchResult")
mapping_file = args.mapping
output_subfolder = args.output


def safe_filename(name):
    name = name.lower()
    name = re.sub(r'\W+', '_', name)
    return name.strip('_')

if not os.path.exists(mapping_file):
    print("❌ 缺少映射文件: {}".format(mapping_file))
    sys.exit(1)

map_df = pd.read_csv(mapping_file, sep='\t', dtype=str).dropna(subset=['organism', 'rename'])

# 构建 rename → organism 映射
rename_to_organism = {}
for i, row in map_df.iterrows():
    for r in str(row['rename']).split(','):
        r = r.strip()
        if r not in rename_to_organism:
            rename_to_organism[r] = set()
        rename_to_organism[r].add(row['organism'].strip())

# 遍历样本文件夹
for sample_folder in os.listdir(tcrmatch_dir):
    sample_path = os.path.join(tcrmatch_dir, sample_folder)
    if not os.path.isdir(sample_path):
        continue

    input_file = os.path.join(sample_path, 'merge_renamed.tsv')
    if not os.path.exists(input_file):
        continue

    try:
        df = pd.read_csv(input_file, sep='\t', dtype=str)
        if df.empty or 'renamed_organism' not in df.columns or 'sequence' not in df.columns:
            continue

        # 加载 trimmed_freq 文件
        trimmed_freq_file = os.path.join(sample_path, 'trimmed_freq.tsv')
        if not os.path.exists(trimmed_freq_file):
            print("⚠️ 缺失 trimmed_freq 文件: {}".format(trimmed_freq_file))
            continue

        trimmed_df = pd.read_csv(trimmed_freq_file, sep='\t', dtype=str)
        if not all(col in trimmed_df.columns for col in ['sequence', 'freq', 'count']):
            print("⚠️ trimmed_freq 文件列不完整: {}".format(trimmed_freq_file))
            continue

        for target_rename in target_renames:
            # 获取 rename 列表
            category_matches = map_df[map_df['category'] == target_rename]
            if not category_matches.empty:
                rename_list = sorted(set(category_matches['rename'].tolist()))
                print("✅ '{}' → 对应 rename: {}".format(target_rename, rename_list))
            else:
                rename_list = [target_rename]
                print("⚠️ 未找到分类，使用默认 rename: {}".format(rename_list))

            # 匹配 sequence 到 organism
            matched = {}
            for i, row in df.iterrows():
                rename_items = [x.strip() for x in str(row['renamed_organism']).split(',')]
                sequence = row['sequence']
                for r in rename_items:
                    if r in rename_list:
                        orgs = rename_to_organism.get(r, [])
                        if sequence not in matched:
                            matched[sequence] = set()
                        matched[sequence].update(orgs)

            if not matched:
                print("ℹ️ 无匹配序列: {} @{}".format(sample_folder, target_rename))
                continue

            # 构建命中 DataFrame
            records = []
            for seq, orgs in matched.items():
                records.append({'sequence': seq, 'organism': ",".join(sorted(orgs))})
            df_hits = pd.DataFrame(records)

            # 合并 freq，不包含 count
            final_df = pd.merge(df_hits, trimmed_df[['sequence', 'freq']], on='sequence', how='left')

            # 对 sequence 去重
            final_df = final_df.drop_duplicates(subset=['sequence'])

            # 输出目录
            output_dir = os.path.join(sample_path, output_subfolder)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # 输出文件路径
            output_file = os.path.join(output_dir, "{}_matched_{}.tsv".format(
                sample_folder, safe_filename(target_rename)))

            # 仅保留需要的列并输出
            final_df[['sequence', 'organism', 'freq']].to_csv(
                output_file, sep='\t', index=False, float_format=None)

            print("✅ 样本完成: {} @{} → {}".format(sample_folder, target_rename, output_file))


    except Exception as e:
        print("❌ 错误 {}: {}".format(sample_folder, e))
