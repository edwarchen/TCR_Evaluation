#!/bin/bash

# Description: 循环后台运行多个 batch 的 TCRmatchpipeline

# 批次列表
batch_list=(
    230903_TCRcell_TB
    231209_TCRcell_TB
    250303_bdszB_TCR
    250529PCRBias_TCR
    250716_TCRcell
    230828_TCRcell_TB
    230921_TCRcell_TB
    250121_TCRcell
    250519_TCRcell
    250712_TCRcell
)

# 日志目录
LOG_DIR="/haplox/users/xuliu/TCR_Project/Results/TCRmatchlogs"
mkdir -p "$LOG_DIR"

# 并发数
MAX_JOBS=3

# 循环批次执行 pipeline
for batch_id in "${batch_list[@]}"; do
    echo "Starting batch: $batch_id"
    nohup sh /haplox/users/xuliu/TCR_Project/scripts/TCRmatchscripts/TCRmatchpipeline.sh "$batch_id" \
        > "$LOG_DIR/${batch_id}.log" 2>&1 &

    # 控制同时运行的任务数
    while (( $(jobs -r | wc -l) >= MAX_JOBS )); do
        sleep 5
    done
done

wait
echo "All batches finished."
