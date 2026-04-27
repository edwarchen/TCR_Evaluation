#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
from tcrdist.repertoire import TCRrep
from Levenshtein import distance
import types
import warnings

# mock statsmodels 避免依赖报错
fake_sm = types.ModuleType("statsmodels")
fake_sm.api = types.ModuleType("statsmodels.api")
sys.modules["statsmodels"] = fake_sm
sys.modules["statsmodels.api"] = fake_sm.api

warnings.filterwarnings("ignore")

# -------------------- 配置 --------------------
ROOT_DIR = "/haplox/users/xuliu/TCR_Project"
OUT_DIR_NAME = "TCRdist_output"
DB_PATH = "/haplox/users/xuliu/miniconda3/envs/tcrdist3/lib/python3.9/site-packages/tcrdist/db/alphabeta_gammadelta_db.tsv"
MAX_SEQUENCES = 5000  # 新增：最大序列数限制

# -------------------- 查找 _formatted.tsv 文件 --------------------
def find_formatted_files(relative_dirs):
    files = []
    for rel_dir in relative_dirs:
        abs_dir = os.path.join(ROOT_DIR, rel_dir)
        if not os.path.exists(abs_dir):
            print(f"[WARN] 目录不存在: {abs_dir}")
            continue
        for root, _, fn_list in os.walk(abs_dir):
            for fn in fn_list:
                if fn.endswith("_formatted.tsv"):
                    files.append(os.path.join(root, fn))
    return files

# -------------------- 读取数据库中的合法 V/J --------------------
def load_valid_genes(db_path):
    db = pd.read_csv(db_path, sep="\t", dtype=str)
    valid_v = set(db.loc[db['region']=='V', 'id'])
    valid_j = set(db.loc[db['region']=='J', 'id'])
    return valid_v, valid_j

VALID_V, VALID_J = load_valid_genes(DB_PATH)

# -------------------- 核心处理 --------------------
def run_tcrdist(path):
    sample = os.path.basename(path).replace("_formatted.tsv", "")
    out_dir = os.path.join(os.path.dirname(path), OUT_DIR_NAME)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n[INFO] ==== 处理样本: {sample} ====")

    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
        df.columns = df.columns.str.strip()  # 去掉可能的空格
    except Exception as e:
        print(f"[WARN] 文件读取失败: {e}")
        return

    # 支持实际列名自动映射
    col_map = {
        "AASeq": "cdr3aa",
        "NNSeq": "cdr3nt",
        "Vregion": "v",
        "Jregion": "j",
        "cloneCount": "count",
        "cloneFraction": "freq"
    }
    df = df.rename(columns=col_map)

    required_cols = ["cdr3aa", "cdr3nt", "v", "j", "count", "freq"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"[WARN] {sample} 缺少必要列: {missing}，跳过")
        return

    # -------------------- 检查合法性 --------------------
    df['valid_cdr3'] = df['cdr3aa'].apply(
        lambda x: pd.notna(x) and len(x) > 0 and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in x.upper())
    )
    df['valid_v'] = df['v'].isin(VALID_V)
    df['valid_j'] = df['j'].isin(VALID_J)

    # 保存被过滤的行
    filtered_out = df[~(df['valid_cdr3'] & df['valid_v'] & df['valid_j'])]
    if not filtered_out.empty:
        filtered_out.to_csv(os.path.join(out_dir, f"{sample}_filtered_out.tsv"), sep="\t", index=False)
        print(f"[INFO] 已保存被过滤行: {sample}_filtered_out.tsv")

    # 保留合法行
    df = df[df['valid_cdr3'] & df['valid_v'] & df['valid_j']]

    if df.empty:
        print(f"[WARN] {sample} 没有合法行，跳过 TCRdist")
        return

    # 强制数值列
    df['count'] = pd.to_numeric(df['count'], errors='coerce').fillna(0).astype(int)
    df['freq'] = pd.to_numeric(df['freq'], errors='coerce').fillna(0.0).astype(float)

    # 新增：按count降序排序并取前MAX_SEQUENCES个序列
    original_count = len(df)
    if original_count > MAX_SEQUENCES:
        df = df.sort_values('count', ascending=False).head(MAX_SEQUENCES)
        print(f"[INFO] 序列数从 {original_count} 限制到 {MAX_SEQUENCES} (按count降序)")
    else:
        print(f"[INFO] 序列数: {original_count} (未超过限制)")

    # 重命名列
    df = df.rename(columns={
        "cdr3aa": "cdr3_b_aa",
        "v": "v_b_gene",
        "j": "j_b_gene"
    })

    df = df.dropna(subset=["cdr3_b_aa","v_b_gene","j_b_gene"])
    df["unique_id"] = df["cdr3_b_aa"] + "_" + df["v_b_gene"].astype(str) + "_" + df["j_b_gene"].astype(str)
    df = df.drop_duplicates("unique_id")
    df = df.set_index("unique_id")

    # 保存 unique_id_map.csv
    df[["cdr3_b_aa","v_b_gene","j_b_gene","freq"]].to_csv(os.path.join(out_dir,f"{sample}_unique_id_map.csv"))
    print(f"[OK] 已生成 unique_id_map.csv")

    # -------------------- 运行 TCRdist --------------------
    try:
        tr = TCRrep(
            cell_df=df,
            organism="human",
            chains=["beta"],
            db_file=DB_PATH,
            compute_distances=True
        )
        
        # 修复：检查可用的距离矩阵属性
        dist_matrix = None
        if hasattr(tr, 'pw_beta'):
            dist_matrix = tr.pw_beta
            print("[INFO] 使用 pw_beta 距离矩阵")
        elif hasattr(tr, 'pws_beta'):
            dist_matrix = tr.pws_beta
            print("[INFO] 使用 pws_beta 距离矩阵")
        elif hasattr(tr, 'rw_beta'):
            dist_matrix = tr.rw_beta
            print("[INFO] 使用 rw_beta 距离矩阵")
        else:
            # 尝试通过调试找到正确的属性名
            print("[DEBUG] 可用的TCRrep属性:")
            for attr in dir(tr):
                if 'beta' in attr.lower() or 'pw' in attr.lower() or 'dist' in attr.lower():
                    print(f"  - {attr}: {type(getattr(tr, attr))}")
            raise AttributeError("未找到可用的距离矩阵属性")
        
        # 创建距离矩阵DataFrame
        dist = pd.DataFrame(dist_matrix, index=df.index, columns=df.index)
        outp = os.path.join(out_dir,f"{sample}_tcrdist_matrix.csv")
        dist.to_csv(outp)
        print(f"[OK] TCRdist matrix 已生成: {outp}")

    except Exception as e:
        print(f"[!!! ALERT !!!] TCRdist 失败: {e}")
        print("[INFO] 调用 Levenshtein fallback，仅基于 CDR3")
        seqs = df["cdr3_b_aa"].tolist()
        ids = df.index.tolist()
        n = len(ids)
        mat = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                d = distance(seqs[i], seqs[j])
                mat[i][j] = mat[j][i] = d
        dist = pd.DataFrame(mat, index=ids, columns=ids)
        outp = os.path.join(out_dir,f"{sample}_cdr3only_levenshtein.csv")
        dist.to_csv(outp)
        print(f"[!!! ALERT !!!] Levenshtein matrix 已生成: {outp}")

# -------------------- 主程序 --------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使用方法：python TCRdist_batch_cleaned.py <相对目录1> [相对目录2 ...]")
        sys.exit(1)

    relative_dirs = sys.argv[1:]
    files = find_formatted_files(relative_dirs)
    print(f"[INFO] 共检测到 {len(files)} 个 formatted 文件")

    for f in files:
        run_tcrdist(f)

    print("\n[DONE] 所有批次处理完成")