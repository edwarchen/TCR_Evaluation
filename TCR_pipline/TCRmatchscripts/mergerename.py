# -*- coding: utf-8 -*-
import os
import pandas as pd
import argparse


def rename_allmerge_organism(base_dir, mapping_file, organism_col='organism'):
    # 读取映射表
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    mapping_df['organism'] = mapping_df['organism'].astype(str).str.strip()
    mapping_df['rename'] = mapping_df['rename'].astype(str).str.strip()

    organism_to_rename = dict(
        zip(mapping_df['organism'], mapping_df['rename'])
    )

    def map_organisms(org_str):
        if pd.isna(org_str) or str(org_str).strip() == '':
            return ''

        organisms = [o.strip() for o in str(org_str).split(',')]
        renamed = [organism_to_rename.get(o, o) for o in organisms]
        return ','.join(renamed)

    for sample_dir in os.listdir(base_dir):
        sample_path = os.path.join(base_dir, sample_dir)
        if not os.path.isdir(sample_path):
            continue

        allmerge_file = os.path.join(sample_path, 'allmerge.tsv')
        if not os.path.exists(allmerge_file):
            continue

        try:
            df = pd.read_csv(allmerge_file, sep='\t')

            if organism_col not in df.columns:
                print(f"❌ {allmerge_file} 中未找到列: {organism_col}")
                continue

            df['renamed_organism'] = df[organism_col].apply(map_organisms)

            output_file = os.path.join(sample_path, 'merge_renamed.tsv')
            df.to_csv(output_file, sep='\t', index=False)
            print(f"✅ Processed: {output_file}")

        except Exception as e:
            print(f"❌ Error processing {allmerge_file}: {e}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', required=True)
    parser.add_argument('--mapping_file', required=True)
    parser.add_argument(
        '--organism_col',
        default='organism',
        help='organism 列名（默认 organism）'
    )
    args = parser.parse_args()

    rename_allmerge_organism(
        args.base_dir,
        args.mapping_file,
        args.organism_col
    )
