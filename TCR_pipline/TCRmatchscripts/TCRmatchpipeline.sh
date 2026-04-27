#!/bin/bash

set -e  # 一旦出错就退出脚本
set -o pipefail

#脚本使用：bash TCRmatchpipeline.sh 250712_TCRcell <output file>

# 检查是否提供了项目名
if [ $# -lt 1 ]; then
    echo "用法: $0 <ProjectName>"
    echo "例如: $0 250712_TCRcell"
    exit 1
fi

PROJECT=$1

# 基础路径
BASE_DIR="/haplox/users/xuliu/TCR_Project/Results/${PROJECT}"
TCRscripts='/haplox/users/xuliu/TCR_Project/scripts/TCRmatchscripts'
TCRMatch='/haplox/users/xuliu/TCR_Project/scripts/TCRMatch-1.3/tcrmatch'
TCRMatchDB='/haplox/users/xuliu/TCR_Project/scripts/TCRMatch-1.3/data/trimmed_cleaned_anti.tsv'

# 配置路径
TCRmatchpath="${BASE_DIR}/TCRmatchResult"
metadata="${BASE_DIR}/metadata.tsv"
name_file="${TCRscripts}/combined_name.tsv"
output_dir="${BASE_DIR}"
cutoff="0.97"


echo "== Step 1: TCRmatch =="
bash /haplox/users/xuliu/TCR_Project/scripts/TCRmatchscripts/TCRmatch.sh $metadata $output_dir $TCRMatch $cutoff $TCRMatchDB  #cutoff值

echo "== Step 2: Processing & Rename =="
python $TCRscripts/mergefreq.py --input_dir $TCRmatchpath
python $TCRscripts/mergerename.py --base_dir $TCRmatchpath --mapping_file $name_file
python $TCRscripts/splitorganism.py --base_dir $TCRmatchpath
python $TCRscripts/countrename.py --base_dir $TCRmatchpath

echo "== Step 3: Combine_all_sample =="
python $TCRscripts/combine_all_counts.py --base_dir $TCRmatchpath

echo "== ✅ Pipeline Completed Successfully =="