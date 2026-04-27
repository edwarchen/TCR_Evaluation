# -*- coding: utf-8 -*-
import os
import csv
from collections import defaultdict, OrderedDict
from decimal import Decimal
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--base_dir', required=True, help='Base directory for TCR results')
args = parser.parse_args()
base_dir = args.base_dir

output_file = os.path.join(base_dir, "final_combined.tsv")

all_organisms = set()
sample_data = {}

# 收集所有 rename，读取每个样本的 rename_summary.tsv
for sample_folder in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample_folder)
    if not os.path.isdir(sample_path):
        continue

    input_file = os.path.join(sample_path, "organism_summary.tsv")
    if not os.path.exists(input_file):
        continue

    sample_name = sample_folder
    data = {}
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            organism = row['organism']
            freq = Decimal(row['freq'])
            count = int(row['count'])
            data[organism] = {'freq': freq, 'count': count}
            all_organisms.add(organism)

    sample_data[sample_name] = data

# 构建列名：organism_freq, organism_count
sorted_organisms = sorted(list(all_organisms))
columns = []
for r in sorted_organisms:
    columns.append(r + "_freq")
    columns.append(r + "_count")

# 写入 final_combined.tsv
with open(output_file, 'w') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["sample"] + columns)

    for sample in sorted(sample_data.keys()):
        row = [sample]
        for r in sorted_organisms:
            entry = sample_data[sample].get(r, {'freq': Decimal('0.0'), 'count': 0})
            row.append("{:.15e}".format(entry['freq']))  # 科学计数法
            row.append(str(entry['count']))
        writer.writerow(row)

print("✅ Completed: {}".format(output_file))
