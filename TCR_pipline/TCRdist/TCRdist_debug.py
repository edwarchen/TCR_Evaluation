#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import types
import warnings
import pandas as pd
from tcrdist.repertoire import TCRrep

# =================================================
# 解决 tcrdist 的 statsmodels 旧依赖问题
# =================================================
fake_sm = types.ModuleType("statsmodels")
fake_sm.api = types.ModuleType("statsmodels.api")
sys.modules["statsmodels"] = fake_sm
sys.modules["statsmodels.api"] = fake_sm.api

warnings.filterwarnings("ignore")

# =================================================
# 常量配置
# =================================================
DB_PATH = (
    "/haplox/users/xuliu/miniconda3/envs/tcrdist3/"
    "lib/python3.9/site-packages/tcrdist/db/"
    "alphabeta_gammadelta_db.tsv"
)

OUT_DIR_NAME = "TCRdist_outputtest"
AA_SET = set("ACDEFGHIKLMNPQRSTVWY")

# =================================================
# 工具函数
# =================================================
def sample_name_from_path(path):
    base = os.path.basename(path)
    return base.replace("_formatted.tsv", "").replace("convert.", "")


def load_valid_vj(db_path):
    db = pd.read_csv(db_path, sep="\t", dtype=str)
    valid_v = set(db.loc[db["region"] == "V", "id"])
    valid_j = set(db.loc[db["region"] == "J", "id"])
    return valid_v, valid_j


def is_valid_cdr3aa(seq):
    if not isinstance(seq, str) or len(seq) == 0:
        return False
    return all(c in AA_SET for c in seq)


def is_valid_nt(seq):
    if not isinstance(seq, str) or len(seq) == 0:
        return False
    return set(seq.upper()) <= set("ACGT")


# =================================================
# 主流程
# =================================================
def run_single_sample(formatted_tsv):

    if not os.path.exists(formatted_tsv):
        raise FileNotFoundError(formatted_tsv)

    sample = sample_name_from_path(formatted_tsv)
    outdir = os.path.join(os.path.dirname(formatted_tsv), OUT_DIR_NAME)
    os.makedirs(outdir, exist_ok=True)

    print(f"\n=== 样本：{sample} ===")

    # -----------------------------
    # 读取数据
    # -----------------------------
    df = pd.read_csv(formatted_tsv, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    # 列名统一
    col_map = {
        "cdr3AA": "cdr3aa",
        "AASeq": "cdr3aa",
        "cdr3NT": "cdr3nt",
        "NNSeq": "cdr3nt",
        "V": "v",
        "Vregion": "v",
        "J": "j",
        "Jregion": "j",
        "cloneCount": "count",
        "cloneFraction": "freq"
    }
    df = df.rename(columns=col_map)

    required = ["cdr3aa", "cdr3nt", "v", "j"]
    for c in required:
        if c not in df.columns:
            raise RuntimeError(f"缺少必要列: {c}")

    # -----------------------------
    # 规范化过滤
    # -----------------------------
    VALID_V, VALID_J = load_valid_vj(DB_PATH)

    df["valid_cdr3aa"] = df["cdr3aa"].apply(is_valid_cdr3aa)
    df["valid_cdr3nt"] = df["cdr3nt"].apply(is_valid_nt)
    df["valid_v"] = df["v"].isin(VALID_V)
    df["valid_j"] = df["j"].isin(VALID_J)

    filtered_out = df[
        ~(df["valid_cdr3aa"] & df["valid_cdr3nt"] & df["valid_v"] & df["valid_j"])
    ]

    if not filtered_out.empty:
        filtered_out.to_csv(
            os.path.join(outdir, f"{sample}_filtered_out.tsv"),
            sep="\t",
            index=False
        )
        print(f"[INFO] 过滤掉 {len(filtered_out)} 条非规范序列")

    df = df[
        df["valid_cdr3aa"] & df["valid_cdr3nt"] & df["valid_v"] & df["valid_j"]
    ].copy()

    if df.empty:
        raise RuntimeError("过滤后无任何合法序列，终止")

    # -----------------------------
    # 转为 TCRrep 所需格式
    # -----------------------------
    df = df.rename(columns={
        "cdr3aa": "cdr3_b_aa",
        "cdr3nt": "cdr3_b_nt",
        "v": "v_b_gene",
        "j": "j_b_gene"
    })
    # ---------- 关键修复 ----------
    df["count"] = pd.to_numeric(df.get("count", 1), errors="coerce").fillna(1).astype(int)
    if "freq" in df.columns:
        df["freq"] = pd.to_numeric(df["freq"], errors="coerce").fillna(0.0)
# -----------------------------

    df["unique_id"] = (
        df["cdr3_b_aa"] + "_" +
        df["v_b_gene"] + "_" +
        df["j_b_gene"]
    )

    df = df.drop_duplicates("unique_id").set_index("unique_id")

    # 保存预处理结果
    df.to_csv(
        os.path.join(outdir, f"{sample}_preprocessed.tsv"),
        sep="\t"
    )

    print(f"[INFO] 合法 TCR 数量: {len(df)}")

    # -----------------------------
    # 计算 TCRdist（严格模式）
    # -----------------------------
    print("[INFO] 开始计算 TCRdist（beta，全模型，无 fallback）")

    tr = TCRrep(
        cell_df=df,
        organism="human",
        chains=["beta"],
        db_file=DB_PATH,
        compute_distances=True
    )

    if not hasattr(tr, "pw_beta"):
        raise RuntimeError("TCRdist 未生成 pw_beta 距离矩阵")

    dist_df = pd.DataFrame(
        tr.pw_beta,
        index=df.index,
        columns=df.index
    )

    out_path = os.path.join(outdir, f"{sample}_pairwise_tcrdist.csv")
    dist_df.to_csv(out_path)

    print(f"[OK] TCRdist 计算完成 -> {out_path}")


# =================================================
# CLI
# =================================================
if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage:")
        print("  python TCRdist_single_sample_strict_full.py sample_formatted.tsv")
        sys.exit(1)

    run_single_sample(sys.argv[1])
