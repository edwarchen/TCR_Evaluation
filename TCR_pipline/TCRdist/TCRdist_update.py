#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import pandas as pd
from tcrdist.repertoire import TCRrep
from Levenshtein import distance
import types
import warnings
import numpy as np

"""
对固定路径/haplox/users/xuliu/TCR_Project/Results或者其他路径下的clonetype文件计算TCRdist，在此之前需要使用TCRseq_to_tcrdistinput.py规范化clonetype文件得到*formatted.tsv
脚本读取的是路径下的*formatted.tsv文件

使用示例：
python TCR_Project/scripts/TCRdist/TCRdist_update.py 250712_TCRcell
python TCR_Project/scripts/TCRdist/TCRdist_update.py /haplox/users/xuliu/TCR_Project/Results/250712_TCRcell

"""

# mock statsmodels 避免依赖报错（tcrdist 内部 import）
fake_sm = types.ModuleType("statsmodels")
fake_sm.api = types.ModuleType("statsmodels.api")
sys.modules["statsmodels"] = fake_sm
sys.modules["statsmodels.api"] = fake_sm.api

warnings.filterwarnings("ignore")

# -------------------- 可配置常量（按需修改） --------------------
ROOT_DIR = "/haplox/users/xuliu/TCR_Project/Results"  # 当你传入的是批次名时，会被拼接到这
OUT_DIR_NAME = "TCRdist_output"
DB_PATH = "/haplox/users/xuliu/miniconda3/envs/tcrdist3/lib/python3.9/site-packages/tcrdist/db/alphabeta_gammadelta_db.tsv"
MAX_SEQUENCES = 5000  # topN 截取（按 count）

# ---------------------------------------------------------------

def find_formatted_files(batch_paths):
    """接受一个批次路径列表（可以是相对批次名，也可以是完整路径），返回所有 *_formatted.tsv 文件（递归）"""
    files = []
    for p in batch_paths:
        # 如果传入的是绝对路径并存在，直接用；否则视作相对于 ROOT_DIR 的批次名
        if os.path.isabs(p) and os.path.exists(p):
            base_dir = p
        else:
            base_dir = os.path.join(ROOT_DIR, p)
        if not os.path.exists(base_dir):
            print(f"[WARN] 批次路径不存在：{base_dir}")
            continue
        for root, _, fn_list in os.walk(base_dir):
            for fn in fn_list:
                if fn.endswith("_formatted.tsv"):
                    files.append(os.path.join(root, fn))
    return sorted(files)


def load_valid_genes(db_path):
    """从 tcrdist db 读取合法 V/J 列表"""
    db = pd.read_csv(db_path, sep="\t", dtype=str)
    valid_v = set(db.loc[db['region']=='V', 'id'])
    valid_j = set(db.loc[db['region']=='J', 'id'])
    return valid_v, valid_j


VALID_V, VALID_J = load_valid_genes(DB_PATH)


def compute_levenshtein_matrix(df, seq_col):
    """基于给定列（氨基酸序列）计算 Levenshtein 距离矩阵，返回 DataFrame"""
    seqs = df[seq_col].tolist()
    ids = list(df.index)
    n = len(seqs)
    mat = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            d = distance(seqs[i], seqs[j])
            mat[i, j] = d
            mat[j, i] = d
    return pd.DataFrame(mat, index=ids, columns=ids)


def safe_sample_name_from_path(path):
    """从文件名里提取更友好的样本名：去掉 convert. 前缀并去掉 _formatted.tsv"""
    base = os.path.basename(path)
    name = base.replace("_formatted.tsv", "")
    name = name.replace("convert.", "")
    return name


def run_single_formatted(path, fallback_samples):
    """处理单个 formatted.tsv：过滤 -> topN -> TCRrep try -> fallback if needed -> 输出矩阵与其他文件"""
    sample = safe_sample_name_from_path(path)
    sample_outdir = os.path.join(os.path.dirname(path), OUT_DIR_NAME)
    os.makedirs(sample_outdir, exist_ok=True)

    print(f"\n--- 处理样本: {sample}")

    # 读入并尝试列名映射
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
        df.columns = df.columns.str.strip()
    except Exception as e:
        print(f"[ERROR] 无法读取 {path}: {e}")
        return

    # 支持常见列名映射（如果你的列名就是下面这些，映射会自动触发）
    col_map = {
        "AASeq": "cdr3aa", "cdr3AA": "cdr3aa", "cdr3_aa": "cdr3aa",
        "NNSeq": "cdr3nt", "cdr3_nt": "cdr3nt", 
        "Vregion": "v", "v_call": "v", "V": "v",
        "Jregion": "j", "j_call": "j", "cloneCount": "count", "cloneFraction": "freq"
    }
    df = df.rename(columns=col_map)

    # 必需列检查
    required = ["cdr3aa", "cdr3nt", "v", "j", "count", "freq"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"[WARN] {sample} 缺少必要列 {missing}，跳过该样本")
        return

    # 过滤合法性
    df['valid_cdr3'] = df['cdr3aa'].apply(
        lambda x: pd.notna(x) and len(x) > 0 and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in x.upper())
    )
    df['valid_v'] = df['v'].isin(VALID_V)
    df['valid_j'] = df['j'].isin(VALID_J)

    filtered_out = df[~(df['valid_cdr3'] & df['valid_v'] & df['valid_j'])]
    if not filtered_out.empty:
        fo_path = os.path.join(sample_outdir, f"{sample}_filtered_out.tsv")
        filtered_out.to_csv(fo_path, sep="\t", index=False)
        print(f"[INFO] {sample} 保存被过滤行 -> {os.path.basename(fo_path)} (共 {len(filtered_out)} 条)")

    df = df[df['valid_cdr3'] & df['valid_v'] & df['valid_j']]
    if df.empty:
        print(f"[WARN] {sample} 过滤后无合法序列，跳过")
        return

    # 强制数值类型
    df['count'] = pd.to_numeric(df['count'], errors='coerce').fillna(0).astype(int)
    df['freq']  = pd.to_numeric(df['freq'], errors='coerce').fillna(0.0).astype(float)

    # 取 topN （按 count）
    original_n = len(df)
    if original_n > MAX_SEQUENCES:
        df = df.sort_values("count", ascending=False).head(MAX_SEQUENCES).copy()
        print(f"[INFO] {sample} 序列 {original_n} -> top{MAX_SEQUENCES} (按 count)")
    else:
        df = df.sort_values("count", ascending=False).copy()
        print(f"[INFO] {sample} 序列 {original_n} (<={MAX_SEQUENCES}) 全部使用")

    # 规范列名以供 TCRrep 使用
    df = df.rename(columns={"cdr3aa": "cdr3_b_aa", "v": "v_b_gene", "j": "j_b_gene"})
    df = df.dropna(subset=["cdr3_b_aa", "v_b_gene", "j_b_gene"])

    # unique_id
    df["unique_id"] = df["cdr3_b_aa"] + "_" + df["v_b_gene"].astype(str) + "_" + df["j_b_gene"].astype(str)
    df = df.drop_duplicates("unique_id").set_index("unique_id")

    # 保存 preprocessed（可追溯）
    preproc_path = os.path.join(sample_outdir, f"{sample}_preprocessed.tsv")
    df.to_csv(preproc_path, sep="\t")
    print(f"[INFO] {sample} 保存预处理表 -> {os.path.basename(preproc_path)} (rows={len(df)})")

    # 尝试运行 TCRrep（TCRdist）
    try:
        print(f"[INFO] {sample} 尝试运行 TCRrep / 计算 TCRdist ...")
        tr = TCRrep(cell_df=df, organism="human", chains=["beta"], db_file=DB_PATH, compute_distances=True)

        # 选取可用矩阵来源（优先 pw_beta；若有 alpha/both 则按需要组合）
        dist_matrix = None
        if hasattr(tr, 'pw_beta'):
            dist_matrix = tr.pw_beta
        elif hasattr(tr, 'pw'):
            dist_matrix = tr.pw
        elif hasattr(tr, 'pws_beta'):
            dist_matrix = tr.pws_beta
        else:
            # 尝试列出候选属性供 debug
            cand = [a for a in dir(tr) if 'pw' in a or 'dist' in a or 'beta' in a or 'alpha' in a]
            raise RuntimeError(f"找不到 pw 距离矩阵属性，候选属性: {cand}")

        dist_df = pd.DataFrame(dist_matrix, index=df.index, columns=df.index)
        out_dist_path = os.path.join(sample_outdir, f"{sample}_pairwise_tcrdist.csv")
        dist_df.to_csv(out_dist_path)
        print(f"[OK] {sample} 成功写出 TCRdist 矩阵 -> {os.path.basename(out_dist_path)}")

    except Exception as e:
        # 进入 fallback
        print(f"[WARN] {sample} TCRdist 计算失败（将 fallback 为 Levenshtein），错误信息: {e}")
        fallback_samples.append(sample)  # 记录会在主流程打印
        try:
            lev_df = compute_levenshtein_matrix(df, seq_col="cdr3_b_aa")
            out_lev_path = os.path.join(sample_outdir, f"{sample}_pairwise_levdist.csv")
            lev_df.to_csv(out_lev_path)
            print(f"[OK] {sample} 成功写出 Levenshtein 矩阵 -> {os.path.basename(out_lev_path)}")
        except Exception as e2:
            print(f"[ERROR] {sample} fallback 计算也失败：{e2}")

    # 保存一个简洁的 unique_id_map（用于下游聚类脚本）
    uidmap_path = os.path.join(sample_outdir, f"{sample}_unique_id_map.csv")
    df_info = df.reset_index()[["unique_id", "cdr3_b_aa", "v_b_gene", "j_b_gene", "count", "freq"]]
    df_info.to_csv(uidmap_path, index=False)
    print(f"[INFO] {sample} 保存 unique_id_map -> {os.path.basename(uidmap_path)}")

    return


def main():
    if len(sys.argv) < 2:
        print("Usage: python batch_tcrdist_with_fallback.py <batch_name_or_fullpath> [<batch...>]")
        print("示例：")
        print("  python batch_tcrdist_with_fallback.py 250712_TCRcell")
        print("  或使用绝对路径： python batch_tcrdist_with_fallback.py /haplox/users/xuliu/TCR_Project/Results/250712_TCRcell")
        sys.exit(1)

    batch_args = sys.argv[1:]
    files = find_formatted_files(batch_args)
    if not files:
        print("[ERROR] 未找到任何 _formatted.tsv 文件，检查你输入的批次名或路径是否正确")
        sys.exit(1)

    print(f"[INFO] 找到 {len(files)} 个 formatted 文件，开始逐个处理...")

    # 记录触发 fallback 的样本
    global fallback_samples
    fallback_samples = []

    for f in files:
        try:
            run_single_formatted(f, fallback_samples)
        except Exception as e:
            print(f"[ERROR] 处理 {f} 时发生未捕获异常: {e}")

    # 运行结束，打印 fallback 样本列表（如有）
    print("\n========== 本次运行触发 fallback（Levenshtein）的样本 ==========")
    if len(fallback_samples) == 0:
        print("无，所有样本均计算了 TCRdist ✅")
    else:
        for s in fallback_samples:
            print(" -", s)
    print("=============================================================\n")


if __name__ == "__main__":
    main()
