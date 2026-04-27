#!/bin/bash

set -e  # 一旦报错就停止

samples=(
	"230828_TCRcell_TB"
    "230921_TCRcell_TB"
    "231209_TCRcell_TB"
    "250303_bdszB_TCR"
	"250519_TCRcell"
    "250529PCRBias_TCR"
	"250712_TCRcell"
    "250716_TCRcell"
    "250817_TCRcell"
    "250716_TCRcell_supple"
)

base_path="/haplox/users/xuliu/TCR_Project/Results"
log_file="run.log"

# 清空之前的日志
#> "$log_file"

#echo "==== Running TCRseq_to_tcrdistinput.py for all samples ====" >> "$log_file"
#for sample in "${samples[@]}"; do
#    echo "Processing $sample ..." | tee -a "$log_file"
#    python TCR_Project/scripts/TCRdist/TCRseq_to_tcrdistinput.py "$sample" >> "$log_file" 2>&1
#done

#echo "==== Running TCRdist_update.py for all samples ====" >> "$log_file"
#for sample in "${samples[@]}"; do
#    echo "Processing $sample ..." | tee -a "$log_file"
#    python TCR_Project/scripts/TCRdist/TCRdist_update.py "$sample" >> "$log_file" 2>&1
#done

#echo "==== Running TCRdist_radius_auto.py for all samples ====" >> "$log_file"
#for sample in "${samples[@]}"; do
#    echo "Processing $sample ..." | tee -a "$log_file"
#    python TCR_Project/scripts/TCRdist/TCRdist_radius_auto.py "$base_path/$sample" >> "$log_file" 2>&1
#done

#echo "==== Running TCRdist_select_lines_auto.py for all samples ====" >> "$log_file"
#for sample in "${samples[@]}"; do
#    echo "Processing $sample ..." | tee -a "$log_file"
#    python TCR_Project/scripts/TCRdist/TCRdist_select_lines_auto.py "$base_path/$sample" top10_cluster_freq_sum >> "$log_file" 2>&1
#done

echo "==== TCRdist_radius_TCR_CI.py for all samples ====" >> "$log_file"
for sample in "${samples[@]}"; do
    echo "Processing $sample ..." | tee -a "$log_file"
    python TCR_Project/scripts/TCRdist/TCRdist_radius_TCR_CI.py --batch "$sample" --radius 30 >> "$log_file" 2>&1
done

echo "All tasks done. Logs saved in $log_file."
