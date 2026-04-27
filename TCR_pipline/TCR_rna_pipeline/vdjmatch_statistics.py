#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import sys
import numpy as np
import math
import subprocess

# ===== 参数 =====
JAR = "/x03_haplox/users/donglf/bin/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar"

metadata_file = os.path.abspath(sys.argv[1])
out_dir = os.path.abspath(sys.argv[2])
final_output = os.path.abspath(sys.argv[3])

os.makedirs(out_dir, exist_ok=True)

# ===== 统计函数 =====
def get_freq_stat(lst):
    total, length = np.sum(lst), len(lst)
    if length == 0:
        return 0, 0, 0
    return total, length, total * math.log(length)


# ===== 读取 metadata =====
meta = pd.read_csv(metadata_file, sep="\t")

print(f">>> Total samples: {len(meta)}")

all_results = []

# ===== 主循环 =====
for _, row in meta.iterrows():
    file_path = row["file_name"]
    sample = row["sample_id"]
    chain = row["chain"]

    print(f"\n>>> Processing: {sample} ({chain})")

    # ===== 1. 复制 clonotype 文件到 out_dir（避免路径问题）=====
    base_input = f"{sample}.input.txt"
    new_file = os.path.join(out_dir, base_input)

    if not os.path.exists(new_file):
        os.system(f"cp {file_path} {new_file}")

    # ===== 2. 构造 metadata =====
    tmp_meta = os.path.join(out_dir, f"meta_{sample}.tsv")
    with open(tmp_meta, "w") as f:
        f.write("file_name\tsample_id\n")
        f.write(f"{base_input}\t{sample}\n")

    # ===== 3. 运行 vdjmatch（在 out_dir）=====
    cmd = f"""
    cd {out_dir} && \
    java -Xmx4G -jar {JAR} match \
        -S human \
        -R {chain} \
        --min-epi-size 2 \
        -m {os.path.basename(tmp_meta)} \
        .
    """

    ret = subprocess.run(cmd, shell=True)

    if ret.returncode != 0:
        print(f"⚠️ vdjmatch failed for {sample}")
        continue

    # ===== 4. 定位输出文件 =====
    result_file = os.path.join(out_dir, f"{sample}.txt")

    if not os.path.exists(result_file):
        print(f"⚠️ result not found: {result_file}")
        continue

    # ===== 5. 读取 vdjmatch 结果 =====
    try:
        df = pd.read_csv(result_file, sep="\t")
    except Exception as e:
        print(f"⚠️ failed to read {result_file}: {e}")
        continue

    # ===== 自动识别列 =====
    if "antigen.species" in df.columns:
        antigen_col = "antigen.species"
    elif "species" in df.columns:
        antigen_col = "species"
    else:
        print(f"⚠️ no antigen column in {sample}")
        continue

    if "freq" not in df.columns:
        print(f"⚠️ no freq column in {sample}")
        continue

    df = df[["freq", antigen_col]]
    df.rename(columns={antigen_col: "antigen.species"}, inplace=True)

    # ===== 6. 聚合统计 =====
    grouped = df.groupby("antigen.species").agg({"freq": list}).reset_index()

    grouped["total_freq"], grouped["count"], grouped["score"] = zip(
        *grouped["freq"].map(get_freq_stat)
    )

    # ===== 7. 拆分 sample 和 chain =====
    if sample.endswith("_TRA"):
        real_sample = sample.replace("_TRA", "")
        real_chain = "TRA"
    elif sample.endswith("_TRB"):
        real_sample = sample.replace("_TRB", "")
        real_chain = "TRB"
    else:
        real_sample = sample
        real_chain = chain

    grouped["sample"] = real_sample
    grouped["chain"] = real_chain

    grouped.drop(columns=["freq"], inplace=True)

    all_results.append(grouped)

    # ===== 8. 删除临时 metadata（可选）=====
    os.remove(tmp_meta)


# ===== 9. 合并输出 =====
if not all_results:
    raise ValueError("❌ No valid vdjmatch results found. 请检查输入文件或路径")

final_df = pd.concat(all_results)
final_df.to_csv(final_output, index=False)

print("\n>>> Pipeline finished!")