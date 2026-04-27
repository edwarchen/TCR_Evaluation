#!/usr/bin/env python3
# coding: utf-8

import os
import sys
import pandas as pd
import subprocess

# ---------------- 配置 ----------------
COS_RAW = "cos://sz-hapseq-1313340341/rawfq"
COS_RESTORE = "cos://sz-hapseq-1313340341/restore/rawfq"
COS_TARGET = "cos://sz-hapdeliver-1313340341/bdsz/LungTCR/tcr_rawfq"
# --------------------------------------

def cos_exists(cos_path):
    """检查 COS 路径下是否有 fastq 文件"""
    cmd = ["coscli", "ls", cos_path]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    for line in result.stdout.splitlines():
        if ".fastq.gz" in line:
            return True
    return False

def get_fastq_files(cos_path):
    """获取路径下所有 fastq 文件，返回完整 COS 路径"""
    cmd = ["coscli", "ls", cos_path]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    files = []
    for line in result.stdout.splitlines():
        parts = line.split("|")
        if len(parts) >= 1:
            path = parts[0].strip()
            if path.endswith(".fastq.gz"):
                if not path.startswith("cos://"):
                    path = "cos://sz-hapseq-1313340341/" + path
                files.append(path)
    return files

def main(sample_info_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    restore_batches = set()
    missing_batches = set()

    transfer_cmd_file = os.path.join(output_dir, "transfer_cmd.sh")
    restore_batch_file = os.path.join(output_dir, "restore_batch_list.txt")
    missing_batch_file = os.path.join(output_dir, "missing_batch_list.txt")

    df = pd.read_csv(sample_info_file)

    with open(transfer_cmd_file, "w") as tf:
        for row in df.itertuples():
            batch, sid = row.Sed_ID, row.Sample_ID

            year_month = '20' + batch[:4]
            base_raw = os.path.join(COS_RAW, year_month, f"{batch}_clinic", sid)
            base_restore = os.path.join(COS_RESTORE, year_month, f"{batch}_clinic", sid)

            if cos_exists(base_raw):
                path_used = base_raw
            elif cos_exists(base_restore):
                path_used = base_restore
                restore_batches.add(batch)
            else:
                missing_batches.add(batch)
                print(f"[MISSING] {sid} ({batch}) not found in either path")
                continue

            fastq_files = get_fastq_files(path_used)
            for fq in fastq_files:
                # 提取文件名
                filename = os.path.basename(fq)
                target_path = os.path.join(COS_TARGET, batch, sid, filename)
                tf.write(f"coscli cp {fq} {target_path}\n")

    # 输出批次文件
    with open(restore_batch_file, "w") as f:
        for b in sorted(restore_batches):
            f.write(b + "\n")
    with open(missing_batch_file, "w") as f:
        for b in sorted(missing_batches):
            f.write(b + "\n")

    print("✅ 已生成文件:")
    print("  ", transfer_cmd_file)
    print("  ", restore_batch_file)
    print("  ", missing_batch_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python3 /haplox/users/xuliu/TCR_Project/scripts/coscli_cp_fq_v2.py /haplox/users/xuliu/TCR_Project/SampleList/TCRfq_meta.csv /haplox/users/xuliu/TCR_Project/SampleList/")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
