#!/bin/bash

###
#第一步得到上机样本文件
Sed=20260304_LH00348_0559_B23JWFJLT4
output=/haplox/users/xuliu/TCR_Project/SampleList/filtered_info.csv
input=/haplox/users/xuliu/TCR_Project/SampleList/info.csv


#scp 从06服务器上下载xlsx文件到36 192.168.1.16:/data/hapreports/hapreports/0.上机信息汇总/上机信息v6.0*share.xlsx

#经过R处理后xlsxtocsv.R得到info.csv
# scp xuliu@192.168.1.36:/home/xuliu/desktop/info.csv TCR_Project/SampleList/ ❗更改传出路径
# 输出表头
echo "Sed_ID,Sample_ID,Lib_number,Name,Panel,Data_ID" > "$output"

# 筛选并输出内容
awk -F',' -v id=$Sed '($1 == id) && ($12 ~ /Adapter_TCR_/) && ($8 ~ /肺结节/) { print $1","$2","$3","$4","$16","$24 }' $input >> $output  ## 
