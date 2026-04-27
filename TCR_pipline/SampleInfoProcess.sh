#!/bin/bash

Sed=250525_A00313_0180_AHTLV5DSXC
output=/haplox/users/xuliu/TCR_Project/SampleList/250525_TCRcell_info.csv
OtherGrep=XKY

echo "Sed_ID,Sample_ID,Lib_number,Name,Panel,Data_ID" > $output
cat /haplox/users/guowei/TCR_Project/scripts/latest_sample.info.tsv | grep -P  $Sed -w | grep  'Adapter_TCR' | grep $OtherGrep |cut -f1,2,3,4,16,24 | sed 's/\t/,/g' >> $output #筛选行和列，对制表符进行全局替换成 ，（即将tsv文件替换成csv文件
