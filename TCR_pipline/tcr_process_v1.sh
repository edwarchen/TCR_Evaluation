#!/bin/bash

batch_id=BCRsample_T3   # Modify according to the analysis content
rawfq_dir_total=/haplox/rawfq/TCR   # default pathway
working_dir=/haplox/users/xuliu/TCR_Project/Results  
sample_info_file=/haplox/users/xuliu/TCR_Project/SampleList/filtered_info.csv #prepare a sample information table
threads=5


rawfq_dir_batch=$rawfq_dir_total/$batch_id
working_dir_batch=$working_dir/$batch_id
mkdir -p $working_dir_batch
mkdir -p $rawfq_dir_batch




function download(){
    rawfq_dir_batch=$1  #路径
    sample_info_file=$2  #样本信息文件

cd $1
/x03_haplox/users/donglf/common_tools/download_scripts/download_fq.py $2 $1 #
bash download_fq.sh #cos目录下载到本地 /haplox/rawfq/TCR/250121_TCRcell 即$1
}



function data_process(){
    sample_info_file=$1  #样本信息文件
    rawfq_dir_batch=$2 #路径
    working_dir_batch=$3 #结果路径 
    threads=$4 #线程

/x03_haplox/users/donglf/haima_pipeline_tools/match_sample.py $sample_info_file $rawfq_dir_batch $rawfq_dir_batch/fq_matched.tsv #样本信息与 FASTQ 文件匹配，生成索引表
cd $working_dir_batch #切换至结果路径
cp $rawfq_dir_batch/fq_matched.tsv . 
cp $sample_info_file . #将匹配文件和样本信息文件复制到当前目录

/x03_haplox/users/donglf/tcr_scripts/tcr_get_shell_fixed_primers.py $working_dir_batch/fq_matched.tsv $working_dir_batch $threads #基于匹配文件生成脚本
bash get_shell.sh | grep sh >run.sh #提取所有包含sh的文件生成run.sh       
/x03_haplox/users/donglf/common_tools/add_wait.py run.sh 8 run2.sh # sample per batch #将run.sh每批8个，分批生成run2
nohup bash run2.sh > run2.log 2>&1 & #后台运行run2
wait #每一步等待后台任务完成
/x03_haplox/users/donglf/tcr_scripts/tcr_qc.py . total_result.tsv  
/x03_haplox/users/donglf/common_tools/fastp_qc.py -i . -o fastp_qc.csv -t .json #将质控文件输出到当前目录
}



#

function data_meta(){
working_dir_batch=$1 #结果路径
cd $working_dir_batch
/x03_haplox/users/donglf/tcr_scripts/get_metadata.py $working_dir_batch metadata.tsv #对应样本和clonetype文件，生成metadata.tsv文件，包含样本信息和每条链的clonetype信息

#random pick 1m reads for statistics, default.
# bash /haplox/users/donglf/tcr_hapyun/stat_from_metadata.sh metadata.tsv $working_dir_batch $working_dir_batch

#random pick 10m reads for statistics. for large sequencing data (>20g)
bash /haplox/users/guowei/TCR_Project/scripts/stat_from_metadata_Samp.sh metadata.tsv $working_dir_batch $working_dir_batch 

}
# samp.sh: tcr_qc.py fastp_qc.py get_total_primer_stat.py get_VJ_stat_from_metadata.py tcr_AA_freq.py tcr_length_freq_group.py tcr_convergence_calculation.py enriched_score_calculation_batch.py cdr3aa_summary.py vdjdb_results_stat.py reshape_vdjdb_stat.py combine_all_stat2.py

###########################################################################
#step1 
download $rawfq_dir_batch $sample_info_file 
#step2
data_process $sample_info_file $rawfq_dir_batch $working_dir_batch $threads
# #step3
data_meta $working_dir_batch  

echo "combined_stat: ${working_dir_batch}/all_combined_stat.csv"
#step4 TCRmatch
sh /haplox/users/xuliu/TCR_Project/scripts/TCRmatchscripts/TCRmatchpipeline.sh $batch_id

#step5 combined all result
python /haplox/users/xuliu/TCR_Project/scripts/TCRmatchscripts/Result_merge.py $working_dir_batch
