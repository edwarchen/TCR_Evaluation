#!/bin/bash

# set -x
metadata=`realpath $1`
output_dir=`realpath $2`
TCRMatch=`realpath $3`
score=$4
TCRMatchDB=`realpath $5`

mkdir -p $output_dir

tail -n +2 "$metadata" | while IFS=$'\t' read -r file_path file_name || [ -n "$file_path" ]; do
  if [ ! -f "$file_path" ]; then
    echo "File not found: $file_path"
    continue
  fi
  sample_dir="$output_dir/TCRmatchResult/${file_name}"
  mkdir -p "$sample_dir"

  # 提取列：sequence（cdr3aa）、trimmed、freq、count
  awk -F'\t' -v OFS='\t' '
  BEGIN {
    cdr_col = 0; freq_col = 0; count_col = 0
    print "sequence", "trimmed", "freq", "count"
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "cdr3aa") cdr_col = i
      if ($i == "freq") freq_col = i
      if ($i == "count") count_col = i
    }
    if (cdr_col == 0 || freq_col == 0 || count_col == 0) {
      print "Error: Required columns not found" > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    seq = $cdr_col
    freq = $freq_col
    count = $count_col
    if (length(seq) >= 5 && seq !~ /[\?\*_]/) {
      trimmed = substr(seq, 2)
      if (substr(trimmed, length(trimmed)) == "F") {
        trimmed = substr(trimmed, 1, length(trimmed)-1)
      }
      print seq, trimmed, freq, count
    }
  }
  ' "$file_path" > "$sample_dir/trimmed_freq_raw.tsv"

  # 聚合同一个 sequence 的 freq 和 count（高精度）
  awk -F'\t' -v OFS='\t' '
  NR == 1 {
    print $1, $2, "freq", "count"
    next
  }
  {
    key = $1 "\t" $2  # sequence + trimmed
    freq[key] += $3
    count[key] += $4
  }
  END {
    for (k in freq) {
      split(k, parts, "\t")
      printf "%s\t%s\t%.16f\t%d\n", parts[1], parts[2], freq[k], count[k]
    }
  }
  ' "$sample_dir/trimmed_freq_raw.tsv" > "$sample_dir/trimmed_freq.tsv"

  # 生成 TCRMatch 输入文件，仅保留 count > 5 的条目
  awk -F'\t' 'NR > 1 && $4 > 5 { print $2 }' "$sample_dir/trimmed_freq.tsv" > "$sample_dir/${file_name}_trimmed_input.txt"
  # 生成 TCRMatch 输入文件
  #cut -f2 "$sample_dir/trimmed_freq.tsv" | tail -n +2 > "$sample_dir/${file_name}_trimmed_input.txt"

  # 运行 TCRMatch
  $TCRMatch -i "$sample_dir/${file_name}_trimmed_input.txt" -t 4 -s $score -d $TCRMatchDB -m 6 > "$sample_dir/${file_name}_tcrmatch.tsv"

  echo "TCRMatch completed: ${file_name}_tcrmatch.tsv"
done

echo "All processing and TCRMatch steps completed."

