#!/bin/bash

#目录下全是.gz文件
SOURCE="/haplox/users/xuliu/TCR_Project/downloadfile2"
#解压后的每个csv整合到一个新的目录
OUTPUT="/haplox/users/xuliu/TCR_Project/downloadfile2/files"

mkdir -p "$OUTPUT"

for f in "$SOURCE"/*.gz; do
    fname=$(basename "$f" .gz)
    gunzip -c "$f" > "$OUTPUT/${fname}"
done
