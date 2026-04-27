#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import re

"""
批次运行示例：
python TCRseq_to_tcrdistinput.py 250716_TCRcell_supple

功能说明：
1. 自动进入 /haplox/users/xuliu/TCR_Project/Results/<批次号>/
2. 自动发现所有 convert.*.clonotypes.TRB.txt 文件
3. 读取 clonetype 文件并进行数据清洗、去重
4. 自动将 v/j 补全等位基因，如 TRBV12-3 → TRBV12-3*01
5. 输出同目录下 *_formatted.tsv
"""

ROOT_DIR = "/haplox/users/xuliu/TCR_Project/Results/"


# -----------------------------
# 修正规因子格式 TRBV12-3 → TRBV12-3*01
# -----------------------------
def fix_allele(gene):
    if pd.isna(gene):
        return gene

    gene = str(gene).strip()

    if "*" in gene:
        return gene   # 已经是等位基因

    m = re.match(r"(TR[AB][VDJ][0-9\-]+)", gene)
    if m:
        return m.group(1) + "*01"

    return gene


# -----------------------------
# 数据准备和去重函数
# -----------------------------
def prepare_and_deduplicate(df, file_path=""):
    """准备数据并去重，返回处理后的DataFrame"""
    if file_path:
        print(f"处理文件: {os.path.basename(file_path)}")
    
    print(f"原始数据: {len(df)} 行, {len(df.columns)} 列")
    
    # 统一列名
    col_map = {}
    for col in df.columns:
        col_lower = str(col).lower()
        if 'cdr3' in col_lower and 'aa' in col_lower:
            col_map[col] = 'cdr3aa'
        elif 'cdr3' in col_lower and 'nt' in col_lower:
            col_map[col] = 'cdr3nt'
        elif col_lower == 'v' or 'vgene' in col_lower or 'vregion' in col_lower:
            col_map[col] = 'v'
        elif col_lower == 'j' or 'jgene' in col_lower or 'jregion' in col_lower:
            col_map[col] = 'j'
        elif 'count' in col_lower or 'clonecount' in col_lower:
            col_map[col] = 'count'
        elif 'freq' in col_lower or 'fraction' in col_lower:
            col_map[col] = 'freq'
    
    if col_map:
        df = df.rename(columns=col_map)
        print(f"重命名后的列: {list(df.columns)}")
    
    # 检查必要列
    required_cols = ['cdr3aa', 'cdr3nt', 'v', 'j']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        print(f"[WARN] 缺少必要列: {missing}，跳过此文件。")
        return None
    
    # 过滤空值
    initial_count = len(df)
    df = df.dropna(subset=['cdr3aa', 'cdr3nt', 'v', 'j'])
    if len(df) < initial_count:
        print(f"移除空值后: {len(df)} / {initial_count}")
    
    # 去除两端空格
    for col in ['cdr3aa', 'cdr3nt', 'v', 'j']:
        df[col] = df[col].astype(str).str.strip()
    
    # 补全 V/J 基因等位基因
    df["v"] = df["v"].apply(fix_allele)
    df["j"] = df["j"].apply(fix_allele)
    
    # 创建唯一标识符并去重
    df['unique_id'] = df['cdr3aa'] + '_' + df['v'] + '_' + df['j']
    
    before_dedup = len(df)
    
    # 对于重复项，合并count和freq
    agg_dict = {
        'cdr3nt': 'first',  # 保留第一个cdr3nt序列
        'count': 'sum',     # count相加
        'freq': 'sum'       # freq相加
    }
    
    # 检查是否有count和freq列，如果没有则添加
    if 'count' not in df.columns:
        df['count'] = 1
    if 'freq' not in df.columns:
        df['freq'] = 0
    
    # 转换为数值类型
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0)
    df["freq"] = pd.to_numeric(df["freq"], errors="coerce").fillna(0)
    
    # 分组聚合去重
    df_unique = df.groupby(['cdr3aa', 'v', 'j'], as_index=False).agg(agg_dict)
    
    after_dedup = len(df_unique)
    
    print(f"去重后唯一序列数: {after_dedup} / {before_dedup} (移除了 {before_dedup - after_dedup} 个重复项)")
    
    # 重新排序列
    output_cols = ["cdr3aa", "cdr3nt", "v", "j", "count", "freq"]
    df_unique = df_unique[output_cols]
    
    # 重置索引
    df_unique = df_unique.reset_index(drop=True)
    df_unique.index.name = 'global_index'
    
    return df_unique


# -----------------------------
# 处理单个 clonotype 文件
# -----------------------------
def process_clonotype_file(infile):
    try:
        # 读取文件
        df = pd.read_csv(infile, sep="\t", dtype=str)
        print(f"[INFO] 读取文件: {os.path.basename(infile)}")
        
        # 使用去重函数处理数据
        df_processed = prepare_and_deduplicate(df, infile)
        
        if df_processed is None:
            return None
        
        return df_processed
        
    except Exception as e:
        print(f"[ERROR] 处理文件 {infile} 时出错: {str(e)}")
        return None


# -----------------------------
# 批量操作
# -----------------------------
def process_batch(batch_name):
    batch_path = os.path.join(ROOT_DIR, batch_name)

    if not os.path.isdir(batch_path):
        print(f"[ERROR] 批次目录不存在：{batch_path}")
        return

    print(f"[INFO] 扫描批次目录：{batch_path}")

    processed_files = 0
    for root, dirs, files in os.walk(batch_path):
        for f in files:
            if f.startswith("convert.") and f.endswith(".clonotypes.TRB.txt"):
                infile = os.path.join(root, f)
                outfile = infile.replace(".txt", "_formatted.tsv")

                df2 = process_clonotype_file(infile)
                if df2 is None:
                    continue

                # 保存处理后的文件
                df2.to_csv(outfile, sep="\t", index=False)
                print(f"[OK] 输出文件：{outfile}")
                processed_files += 1
    
    print(f"[SUMMARY] 共处理 {processed_files} 个文件")


# -----------------------------
# 主程序
# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("使用方法：python TCRseq_to_tcrdistinput.py <批次号>")
        print("示例：python TCRseq_to_tcrdistinput.py 250716_TCRcell_supple")
        sys.exit(1)

    batch_name = sys.argv[1]
    process_batch(batch_name)