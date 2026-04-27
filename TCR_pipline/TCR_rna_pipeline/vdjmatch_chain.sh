#!/bin/bash

set -euo pipefail

JAR="/x03_haplox/users/donglf/bin/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar"
META="$1"
OUT_DIR="$2"

mkdir -p "${OUT_DIR}"
WORK_DIR="${OUT_DIR}/work_tmp"
mkdir -p "${WORK_DIR}"

echo ">>> Start vdjmatch pipeline"

tail -n +2 "$META" | while IFS=$'\t' read -r file sample chain
do
    echo ">>> Processing: ${sample} (${chain})"

    base=$(basename "$file")

    cp "$file" "${WORK_DIR}/${base}"

    TMP_META="meta_${sample}_${chain}.tsv"
    echo -e "file_name\tsample_id" > "${WORK_DIR}/${TMP_META}"
    echo -e "${base}\t${sample}" >> "${WORK_DIR}/${TMP_META}"

    cd "$WORK_DIR"

    out_prefix="${OUT_DIR}/${sample}_${chain}"

    java -Xmx4G -jar "$JAR" match \
        -S human \
        -R ${chain} \
        --min-epi-size 2 \
        -m "$TMP_META" \
        . > "${out_prefix}.log" 2>&1

    result_file=$(ls -t *.txt 2>/dev/null | head -n 1 || true)

    if [ -n "$result_file" ]; then
        mv "$result_file" "${out_prefix}.txt"
        echo "✔ Output: ${out_prefix}.txt"
    else
        echo "⚠️ No result for ${sample} (${chain})"
    fi

    cd - > /dev/null

done

echo ">>> All done!"