#!/bin/bash
set -e

# ==========================
# 参数检查
# ==========================
if [ $# -lt 6 ]; then
echo ""
echo "USAGE:"
echo " $0 input_R1.fastq.gz input_R2.fastq.gz outputdir sample_key threads chain_type"
echo ""
echo " chain_type: TRA / TRB / BOTH"
echo ""
exit 1
fi

R1=$1
R2=$2
outdir=$3
key=$4
thread=$5
CHAIN=$6

shell=$outdir/$key.run.sh

# ==========================
# 基础路径
# ==========================
Bin="/x03_haplox/users/donglf/TCR_chenyr/shell"
Bin2="/x03_haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/"
V_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/V_primer.txt"
J_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/J_primer.txt"
J_rc_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/J_rc_primer.txt"

# ==========================
# 创建目录
# ==========================
mkdir -p $outdir/{Cleanfq,Merge_PE,Map_Clone_Analysis,Stat_Picture,log}

# ==========================
# chain控制
# ==========================
CHAINS_TO_RUN=()

if [ "$CHAIN" == "TRB" ]; then
    CHAINS_TO_RUN=("TRB")
elif [ "$CHAIN" == "TRA" ]; then
    CHAINS_TO_RUN=("TRA")
elif [ "$CHAIN" == "BOTH" ]; then
    CHAINS_TO_RUN=("TRB" "TRA")
else
    echo "ERROR: chain_type must be TRA / TRB / BOTH"
    exit 1
fi

# ==========================
# 写shell
# ==========================
echo "#!/bin/bash" > $shell
echo "set -e" >> $shell
echo "#########################################################################" >> $shell

# ==========================
# fastp
# ==========================
echo "$Bin/fastp-0.19.7/fastp --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -i $R1 -I $R2 -o $outdir/Cleanfq/$key.good_R1.fq.gz -O $outdir/Cleanfq/$key.good_R2.fq.gz -t 1 -T 1 -w $thread -j $outdir/Cleanfq/$key.json -h $outdir/Cleanfq/$key.html >$outdir/log/$key.fastp.log 2>&1" >> $shell

# ==========================
# FLASH
# ==========================
#echo "$Bin/FLASH-1.2.11/FLASH-1.2.11-Linux-x86_64/flash $outdir/Cleanfq/$key.good_R1.fq.gz $outdir/Cleanfq/$key.good_R2.fq.gz -d $outdir/Merge_PE -o $key -m 10 -M 150 -p 33 -r 150 -x 0.1 > $outdir/log/$key.flash.log" >> $shell

# ==========================
# primer筛选
# ==========================
#echo "seqkit grep -s -i --degenerate -f $V_primers $outdir/Merge_PE/$key.extendedFrags.fastq -o $outdir/Merge_PE/$key.extendedFrags.containV.fastq" >> $shell
#echo "perl $Bin/stat_primer.pl -i $outdir/Merge_PE/$key.extendedFrags.fastq -v $outdir/Merge_PE/$key.extendedFrags.containV.fastq -o $outdir/Merge_PE/$key.primer.stat.txt" >> $shell

#echo "seqkit grep -s -i --degenerate -f $J_rc_primers $outdir/Merge_PE/$key.extendedFrags.containV.fastq -o $outdir/Merge_PE/$key.extendedFrags.containVJ.fastq" >> $shell
#echo "perl $Bin/stat_primer.pl -i $outdir/Merge_PE/$key.extendedFrags.containV.fastq -v $outdir/Merge_PE/$key.extendedFrags.containVJ.fastq -o $outdir/Merge_PE/$key.primer.stat.txt" >> $shell

# ==========================
# MiXCR
# ==========================
echo "java -Xmx16g -Xms4g -jar $Bin/mixcr-3.0.4/mixcr.jar align -f --species hs --report $outdir/Map_Clone_Analysis/$key.report -p rna-seq $outdir/Cleanfq/$key.good_R1.fq.gz $outdir/Cleanfq/$key.good_R2.fq.gz $outdir/Map_Clone_Analysis/$key.vdjca -t $thread >$outdir/Map_Clone_Analysis/$key.log 2>&1" >> $shell

echo "java -Xmx16g -Xms4g -jar $Bin/mixcr-3.0.4/mixcr.jar assemble -f -OassemblingFeatures=\"[CDR3]\" $outdir/Map_Clone_Analysis/$key.vdjca $outdir/Map_Clone_Analysis/$key.clna -t $thread >>$outdir/Map_Clone_Analysis/$key.log 2>&1" >> $shell

# ==========================
# export + 后续分析（核心循环）
# ==========================
for c in "${CHAINS_TO_RUN[@]}"; do

echo "echo 'Processing $c ...'" >> $shell

# export
echo "java -Xmx16g -Xms4g -jar $Bin/mixcr-3.0.4/mixcr.jar exportClones -f --chains $c $outdir/Map_Clone_Analysis/$key.clna $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt >>$outdir/Map_Clone_Analysis/$key.log 2>&1" >> $shell

# filter
echo "mv $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.raw.txt" >> $shell
echo "perl $Bin/filter_cdr3_length.pl -i $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.raw.txt -o $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt" >> $shell

# stat
echo "perl $Bin/stat_cdr3_length.pl -i $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt -o $outdir/Stat_Picture/$key.clonotypes.${c}.cdr3.stat.txt -p $outdir/Stat_Picture/$key.cdr3.len.${c}.png >$outdir/log/$key.cdr3.len.${c}.log 2>&1" >> $shell

# vdjtools
echo "java -jar $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S MiXcr $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt $outdir/Map_Clone_Analysis/convert >$outdir/log/$key.convert.${c}.log 2>&1" >> $shell

echo "java -jar $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcBasicStats $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt $outdir/Stat_Picture/$key.${c}.0 >$outdir/log/$key.basic.${c}.log 2>&1" >> $shell

echo "perl $Bin/stat_VJ_count.pl -i $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt -o $outdir/Stat_Picture/$key.${c}.VJ.stat.txt -p $outdir/Stat_Picture/$key.${c}.VJ.heatmap.png >$outdir/log/$key.VJ.${c}.log 2>&1" >> $shell

echo "java -jar $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcSpectratype $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt $outdir/Stat_Picture/$key.${c}.1 >$outdir/log/$key.spectra.${c}.log 2>&1" >> $shell

echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancySpectratype $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt  $outdir/Stat_Picture/$key.2 >$outdir/log/$key.PlotFancySpectratype.${c}.log 2>&1 ">>$shell

echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancyVJUsage  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt  $outdir/Stat_Picture/$key.3  >$outdir/log/$key.PlotFancyVJUsage.${c}.log 2>&1">>$shell

echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotSpectratypeV  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt  $outdir/Stat_Picture/$key.4 >$outdir/log/$key.PlotSpectratypeV.${c}.log 2>&1 ">>$shell



done

echo "perl2 $Bin/total_statistic.v3.pl -i  $outdir/ -o  $outdir/Stat_Picture -k $key -c  $CHAIN  >$outdir/log/$key.total_stat.log 2>&1 ">>$shell


# ==========================
# 清理
# ==========================
echo "rm -rf $outdir/Cleanfq/*fq.gz $outdir/Merge_PE/*fastq" >> $shell

# ==========================
# 提示运行
# ==========================
echo ""
echo "Run with:"
echo "sh $shell > $shell.o 2>&1 &"
echo ""