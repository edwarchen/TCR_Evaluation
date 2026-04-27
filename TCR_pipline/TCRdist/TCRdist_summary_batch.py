#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd


def shannon_entropy(p):
    p = p[p > 0]
    return -np.sum(p * np.log(p))


def summarize_sample(pre_file):
    df = pd.read_csv(pre_file, sep="\t")

    n_seq = df.shape[0]

    # 选择 freq / count
    if "freq" in df.columns:
        freq = df["freq"].astype(float)
    elif "freq_norm" in df.columns:
        freq = df["freq_norm"].astype(float)
    elif "count" in df.columns:
        freq = df["count"].astype(float)
        freq = freq / freq.sum()
    else:
        raise ValueError(f"No freq or count column in {pre_file}")

    freq = freq / freq.sum()

    max_freq = freq.max()
    top10_freq_sum = freq.sort_values(ascending=False).head(10).sum()
    entropy = shannon_entropy(freq.values)

    return {
        "n_seq": n_seq,
        "max_freq": max_freq,
        "top10_freq_sum": top10_freq_sum,
        "entropy": entropy
    }


def summarize_batch(batch_dir):
    records = []

    samples = [
        d for d in os.listdir(batch_dir)
        if os.path.isdir(os.path.join(batch_dir, d))
    ]

    for sample in samples:
        out_dir = os.path.join(
            batch_dir,
            sample,
            "Map_Clone_Analysis",
            "TCRdist_output"
        )

        if not os.path.exists(out_dir):
            continue

        pre_files = [
            f for f in os.listdir(out_dir)
            if f.endswith(".clonotypes.TRB_preprocessed.tsv")
        ]

        if len(pre_files) == 0:
            continue

        pre_file = os.path.join(out_dir, pre_files[0])

        try:
            stats = summarize_sample(pre_file)
        except Exception as e:
            print(f"[WARN] {sample}: {e}")
            continue

        stats["sample"] = sample
        records.append(stats)

    return pd.DataFrame(records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--batch_dir",
        required=True,
        help="Path to batch directory"
    )
    parser.add_argument(
        "--out_file",
        default="batch_sample_summary.tsv"
    )
    args = parser.parse_args()

    df = summarize_batch(args.batch_dir)
    df = df.sort_values("n_seq", ascending=False)

    df.to_csv(args.out_file, sep="\t", index=False)
    print(f"[DONE] Saved summary to {args.out_file}")
