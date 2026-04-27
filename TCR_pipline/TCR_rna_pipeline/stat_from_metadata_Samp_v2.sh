set -x
metadata=`realpath $1`
input_dir=`realpath $2`
output_dir=`realpath $3`
mkdir -p $output_dir && cd $output_dir
/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/tcr_qc_allchain.py $input_dir total_result.tsv
/x03_haplox/users/donglf/common_tools/fastp_qc.py -i $input_dir -o fastp_qc.csv -t .json
#/x03_haplox/users/donglf/tcr_scripts/get_total_primer_stat.py $input_dir total_primer_stat.csv #没有对VJ的对应文件处理
/haplox/users/xuliu/BCR_Project/Script/BCR_script/get_VJ_stat_from_metadata.py $metadata total_V.csv total_J.csv score
/haplox/users/xuliu/BCR_Project/Script/BCR_script/bcr_AA_freq.py $metadata aa_stat.csv
/haplox/users/xuliu/BCR_Project/Script/BCR_script/bcr_length_freq_group.py $metadata length_stat.csv
/haplox/users/xuliu/BCR_Project/Script/BCR_script/bcr_convergence_calculation.py $metadata convergence_stat.csv

awk 'NR==1 {print $0} NR>1 {print $1, $2"_"$3, $3}' OFS="\t" $metadata > metadata_2chains.tsv
mkdir -p vdjtools_stat && cd vdjtools_stat
java -jar /x03_haplox/users/donglf/tcr_hapyun/envs/bin/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcDiversityStats -x 10000000  -m $input_dir/metadata_2chains.tsv . > vdjtools_calstat.log 2>&1 &
cd ..
/x03_haplox/users/donglf/tcr_data_analysis/enriched_cdr3aa/enriched_score_calculation_batch.py metadata_2chains.tsv enrich_stat 10 enriched_cdr3aa.csv #数据库是TRB的结果
bash enrich_stat/run.sh > enrich_stat/run.log 2>&1 &
/x03_haplox/users/donglf/tcr_hapyun/cdr3aa_summary.py metadata_2chains.tsv cdr3aa_stat_summary.csv
mkdir -p vdjmatch && cd vdjmatch
# vdjmatch arguments

# java -Xmx4G -jar /x03_haplox/users/donglf/bin/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRB  --min-epi-size 2 -m /haplox/users/xuliu/TCR_Project/TCR_BCR_RD_Proj/TR_6sample/metadata_2chains.tsv . > vdjmatch.log 2>&1
# /x03_haplox/users/donglf/tcr_hapyun/vdjdb_results_stat.py metadata.txt . vdjdb_anno_summary.csv
#上两步的综合结果
/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/vdjmatch_statistics.py $input_dir/metadata_2chains.tsv . vdjdb_anno_summary.csv
cd ..
/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/reshape_vdjdb_stat_v2.py vdjmatch/vdjdb_anno_summary.csv reshaped_vdjdb_anno_summary.csv
wait
/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/combined_all_stat2.py . all_combined_stat.csv

