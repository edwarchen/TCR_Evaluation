R1=$1
R2=$2
outdir=$3
key=$4
thread=$5
shell=$outdir/$key.run.sh

if [ $# -lt 5  ]; then
echo "";
echo "USAGE:";
echo "	$0 input_R1.fastq.gz input_R2.fastq.gz  outputdir sample_key threads";
echo "";
echo "	input_R1.fastq.gz:	input file xx.R1.fastq.gz";
echo "	input_R2.fastq.gz:	input file xx.R2.fastq.gz";
echo "	outputdir:	output dir ";
echo "	sample_key	a key word of input sample";
echo "	threads:	threads of used";
echo "";
exit 1;
fi
#Bin=$(pwd)
# Bin="$( cd "$( dirname "$0"  )" && pwd  )"
Bin="/x03_haplox/users/donglf/TCR_chenyr/shell" # give the absolute path of the bin
V_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/V_primer.txt"
J_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/J_primer.txt"
J_rc_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/J_rc_primer.txt"
if [ !  -d $outdir ];then  mkdir -p $outdir;  fi
if [ !  -d  $outdir/Cleanfq ];then mkdir -p $outdir/Cleanfq; fi
if [ !  -d  $outdir/Merge_PE ];then mkdir -p $outdir/Merge_PE; fi
if [ !  -d  $outdir/Map_Clone_Analysis ];then mkdir -p $outdir/Map_Clone_Analysis; fi
if [ !  -d  $outdir/Stat_Picture ];then mkdir -p $outdir/Stat_Picture; fi
if [ !  -d  $outdir/log ];then mkdir -p $outdir/log; fi
 
echo "#sh /x03_haplox/users/donglf/tcr_scripts/TCR_analysis_pipeline.v2_1.sh $R1 $R2  $outdir $key $thread ">$shell
echo "######################################################################### ">>$shell
#fastp QC
echo "$Bin/fastp-0.19.7/fastp   --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  -i $R1 -I  $R2 -o $outdir/Cleanfq/$key.good_R1.fq.gz  -O  $outdir/Cleanfq/$key.good_R2.fq.gz -t 1 -T 1 -w $thread -j  $outdir/Cleanfq/$key.json -h $outdir/Cleanfq/$key.html >$outdir/log/$key.fastp.log 2>&1 ">>$shell

echo "$Bin/FLASH-1.2.11/FLASH-1.2.11-Linux-x86_64/flash $outdir/Cleanfq/$key.good_R1.fq.gz $outdir/Cleanfq/$key.good_R2.fq.gz -d $outdir/Merge_PE -o $key -m 10 -M 150 -p 33 -r 150 -x 0.1 > $outdir/log/$key.flash.log ">>$shell

# echo "/x03_haplox/users/donglf/miniconda3/envs/tcr_hapyun/bin/seqkit  grep -s -i -j $thread --degenerate -f $V_primers    $outdir/Cleanfq/$key.good_R1.fq.gz -o $outdir/Cleanfq/$key.good_R1.containV.fq.gz ">>$shell
# echo "perl $Bin/stat_primer.pl  -i $outdir/Cleanfq/$key.good_R1.fq.gz  -v $outdir/Cleanfq/$key.good_R1.containV.fq.gz -o $outdir/Merge_PE/$key.primer.stat.txt  " >>$shell
# echo "/x03_haplox/users/donglf/miniconda3/envs/tcr_hapyun/bin/seqkit  grep -s -i -j $thread --degenerate -f $J_primers    $outdir/Cleanfq/$key.good_R2.fq.gz -o $outdir/Cleanfq/$key.good_R2.containJ.fq.gz ">>$shell
# echo "perl $Bin/stat_primer.pl  -i $outdir/Cleanfq/$key.good_R2.fq.gz  -v $outdir/Cleanfq/$key.good_R2.containJ.fq.gz  -o $outdir/Merge_PE/$key.primer.stat.txt  " >>$shell

#echo "/x03_haplox/users/donglf/miniconda3/envs/tcr_hapyun/bin/seqkit  grep -s -i --degenerate -f $V_primers  $outdir/Merge_PE/$key.extendedFrags.fastq -o $outdir/Merge_PE/$key.extendedFrags.containV.fastq ">>$shell

#echo "perl $Bin/stat_primer.pl  -i $outdir/Merge_PE/$key.extendedFrags.fastq  -v $outdir/Merge_PE/$key.extendedFrags.containV.fastq -o $outdir/Merge_PE/$key.primer.stat.txt  " >>$shell

#echo "/x03_haplox/users/donglf/miniconda3/envs/tcr_hapyun/bin/seqkit  grep -s -i --degenerate -f $J_rc_primers   $outdir/Merge_PE/$key.extendedFrags.containV.fastq -o $outdir/Merge_PE/$key.extendedFrags.containVJ.fastq ">>$shell

#echo "perl $Bin/stat_primer.pl  -i $outdir/Merge_PE/$key.extendedFrags.containV.fastq  -v $outdir/Merge_PE/$key.extendedFrags.containVJ.fastq -o $outdir/Merge_PE/$key.primer.stat.txt  " >>$shell

echo "java -Xmx16g -Xms4g -jar   $Bin/mixcr-3.0.4/mixcr.jar  align  --species hs --report $outdir/Map_Clone_Analysis/$key.report -p rna-seq -OvParameters.geneFeatureToAlign=VGeneWithP -OvParameters.parameters.floatingLeftBound=true -OjParameters.parameters.floatingRightBound=true -OcParameters.parameters.floatingRightBound=false  $outdir/Merge_PE/$key.extendedFrags.fastq  $outdir/Map_Clone_Analysis/$key.vdjca -t $thread >$outdir/Map_Clone_Analysis/$key.log 2>&1 ">>$shell
echo "java -Xmx16g -Xms4g  -jar $Bin/mixcr-3.0.4/mixcr.jar  assemble   --write-alignments -OassemblingFeatures="[CDR3]" -OseparateByV=false -OseparateByJ=false -OseparateByC=false $outdir/Map_Clone_Analysis/$key.vdjca $outdir/Map_Clone_Analysis/$key.clna  -t $thread   >>$outdir/Map_Clone_Analysis/$key.log 2>&1 ">>$shell

echo "java -Xmx16g -Xms4g  -jar $Bin/mixcr-3.0.4/mixcr.jar  exportClones    --chains TRB   $outdir/Map_Clone_Analysis/$key.clna  $outdir/Map_Clone_Analysis/$key.clonotypes.TRB.txt  >>$outdir/Map_Clone_Analysis/$key.log 2>&1 ">>$shell
echo "java -Xmx16g -Xms4g  -jar $Bin/mixcr-3.0.4/mixcr.jar  exportClones    --chains TRA   $outdir/Map_Clone_Analysis/$key.clna  $outdir/Map_Clone_Analysis/$key.clonotypes.TRA.txt  >>$outdir/Map_Clone_Analysis/$key.log 2>&1 ">>$shell

###############################add filter CDR3 length #########################################
echo "mv $outdir/Map_Clone_Analysis/$key.clonotypes.TRB.txt  $outdir/Map_Clone_Analysis/$key.clonotypes.TRB.raw.txt  ">>$shell
echo "mv $outdir/Map_Clone_Analysis/$key.clonotypes.TRA.txt  $outdir/Map_Clone_Analysis/$key.clonotypes.TRA.raw.txt  ">>$shell
echo "perl $Bin/filter_cdr3_length.pl  -i  $outdir/Map_Clone_Analysis/$key.clonotypes.TRB.raw.txt -o  $outdir/Map_Clone_Analysis/$key.clonotypes.TRB.txt">>$shell
echo "perl $Bin/filter_cdr3_length.pl  -i  $outdir/Map_Clone_Analysis/$key.clonotypes.TRA.raw.txt -o  $outdir/Map_Clone_Analysis/$key.clonotypes.TRA.txt">>$shell
###############################statistic and draw pictures ######################################

echo "perl $Bin/stat_contig_length.pl -i  $outdir/Merge_PE/$key.hist    -o  $outdir/Merge_PE/$key.hist.2 -p $outdir/Stat_Picture/$key.contig.len.png     >$outdir/log/$key.contig.len.log 2>&1 ">>$shell

echo "perl $Bin/stat_cdr3_length.pl  -i   $outdir/Map_Clone_Analysis/$key.clonotypes.TRB.txt   -o   $outdir/Stat_Picture/$key.clonotypes.TRB.cdr3.stat.txt  -p $outdir/Stat_Picture/$key.cdr3.len.TRB.png     >$outdir/log/$key.cdr3.len.log 2>&1 ">>$shell
echo "perl $Bin/stat_cdr3_length.pl  -i   $outdir/Map_Clone_Analysis/$key.clonotypes.TRA.txt   -o   $outdir/Stat_Picture/$key.clonotypes.TRA.cdr3.stat.txt  -p $outdir/Stat_Picture/$key.cdr3.len.TRA.png     >$outdir/log/$key.cdr3.len.log 2>&1 ">>$shell

echo "java -jar $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S MiXcr $outdir/Map_Clone_Analysis/$key.clonotypes.TRB.txt $outdir/Map_Clone_Analysis/convert  >$outdir/log/$key.convert.log 2>&1 ">>$shell
echo "java -jar $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S MiXcr $outdir/Map_Clone_Analysis/$key.clonotypes.TRA.txt $outdir/Map_Clone_Analysis/convert  >$outdir/log/$key.convert.log 2>&1 ">>$shell

echo "java -jar $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcBasicStats  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRB.txt   $outdir/Stat_Picture/$key.0  >$outdir/log/$key.CalcBasicStats.log 2>&1">>$shell
echo "java -jar $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcBasicStats  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRA.txt   $outdir/Stat_Picture/$key.0  >$outdir/log/$key.CalcBasicStats.log 2>&1">>$shell

echo "perl  $Bin/stat_VJ_count.pl -i $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRB.txt  -o $outdir/Stat_Picture/$key.TRB.VJ.stat.txt  -p $outdir/Stat_Picture/$key.TRB.VJ.heatmap.png >$outdir/log/$key.stat_VJ_count.log 2>&1" >>$shell
echo "perl  $Bin/stat_VJ_count.pl -i $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRA.txt  -o $outdir/Stat_Picture/$key.TRA.VJ.stat.txt  -p $outdir/Stat_Picture/$key.TRA.VJ.heatmap.png >$outdir/log/$key.stat_VJ_count.log 2>&1" >>$shell

echo "perl  $Bin/draw_3D_VJ_pic.pl   -i $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRB.txt  -o $outdir/Stat_Picture   -k  $key   >$outdir/log/$key.3D.VJ.log 2>&1 " >>$shell
echo "perl  $Bin/draw_3D_VJ_pic.pl   -i $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRA.txt  -o $outdir/Stat_Picture   -k  $key   >$outdir/log/$key.3D.VJ.log 2>&1 " >>$shell

echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcSpectratype  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRB.txt    $outdir/Stat_Picture/$key.TRB.1  >$outdir/log/$key.CalcSpectratype.log 2>&1">>$shell
echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcSpectratype  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRA.txt    $outdir/Stat_Picture/$key.TRA.1  >$outdir/log/$key.CalcSpectratype.log 2>&1">>$shell

echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancySpectratype $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRB.txt  $outdir/Stat_Picture/$key.TRB.2 >$outdir/log/$key.PlotFancySpectratype.log 2>&1 ">>$shell
echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancySpectratype $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRA.txt  $outdir/Stat_Picture/$key.TRA.2 >$outdir/log/$key.PlotFancySpectratype.log 2>&1 ">>$shell

echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancyVJUsage  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRB.txt  $outdir/Stat_Picture/$key.TRB.3  >$outdir/log/$key.PlotFancyVJUsage.log 2>&1">>$shell
echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancyVJUsage  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRA.txt  $outdir/Stat_Picture/$key.TRA.3  >$outdir/log/$key.PlotFancyVJUsage.log 2>&1">>$shell

echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotSpectratypeV  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRB.txt  $outdir/Stat_Picture/$key.TRB.4 >$outdir/log/$key.PlotSpectratypeV.log 2>&1 ">>$shell
echo "java -jar   $Bin/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotSpectratypeV  $outdir/Map_Clone_Analysis/convert.$key.clonotypes.TRA.txt  $outdir/Stat_Picture/$key.TRA.4 >$outdir/log/$key.PlotSpectratypeV.log 2>&1 ">>$shell

echo "perl $Bin/total_statistic.v1.pl -i  $outdir/   -o  $outdir/Stat_Picture -k $key   >$outdir/log/$key.total_stat.log 2>&1 ">>$shell

echo "rm -fr $outdir/haplox  vj_pairing_plot.r  fancy_spectratype.r ">>$shell
echo "rm -fr $outdir/Cleanfq/*fq.gz  $outdir/Merge_PE/*fastq  ">>$shell


echo " "
echo "sh $outdir/$key.run.sh > $outdir/$key.run.sh.o 2>&1 &"
echo " "
