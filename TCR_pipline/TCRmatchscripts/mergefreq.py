import os
import pandas as pd
import glob
import argparse

# Argument parser for command-line execution
parser = argparse.ArgumentParser(description="Merge *.tcrmatch.tsv and trimmed_freq files for each sample folder.")
parser.add_argument(
    '--input_dir',
    required=True,
    help='Path to TCRmatchResult directory (e.g., ../data/TCRmatchResult)'
)
args = parser.parse_args()

base_dir = os.path.abspath(args.input_dir)

for sample_dir in os.listdir(base_dir):
    sample_path = os.path.join(base_dir, sample_dir)
    if os.path.isdir(sample_path):
        # *.tcrmatch.tsv file (only first match is used)
        tsv_files = glob.glob(os.path.join(sample_path, '*tcrmatch.tsv'))
        if not tsv_files:
            print("[WARNING] No .tcrmatch.tsv file found in folder: " + sample_dir)
            continue
        tsv_file = tsv_files[0]

        # trimmed_freq file
        trimmed_freq_path = os.path.join(sample_path, 'trimmed_freq.tsv')
        if not os.path.exists(trimmed_freq_path):
            print("[WARNING] No trimmed_freq file found in folder: " + sample_dir)
            continue

        try:
            tcrmatch_df = pd.read_csv(tsv_file, sep='\t')
            trimmed_df = pd.read_csv(trimmed_freq_path, sep='\t')
        except Exception as e:
            print("[ERROR] Failed to read files in folder: " + sample_dir)
            print("Reason: " + str(e))
            continue

        # Merge on input_sequence (left) and trimmed (right)
        merged_df = pd.merge(tcrmatch_df, trimmed_df, left_on='trimmed_input_sequence', right_on='trimmed', how='inner')
        
        # Delete
        if 'trimmed' in merged_df.columns:
            merged_df = merged_df.drop(columns=['trimmed'])

        # Remove completely duplicate rows
        merged_df = merged_df.drop_duplicates()

        # Identify unmatched rows
        unmatched_df = tcrmatch_df[~tcrmatch_df['trimmed_input_sequence'].isin(merged_df['trimmed_input_sequence'])]
        
        # Remove completely duplicate rows in unmatched_df (optional)
        unmatched_df = unmatched_df.drop_duplicates()

        # Write to file: Contains complete unmatched lines
        unmatched_df.to_csv(os.path.join(sample_path, 'unmatchline.tsv'), sep='\t', index=False)

        # de-duplicated
        unmatched_sequences = unmatched_df['trimmed_input_sequence'].drop_duplicates()

        #print("[INFO] Unique unmatched sequences in sample '{}':".format(sample_dir))
        for seq in unmatched_sequences:
            print(seq)


        merged_df.to_csv(os.path.join(sample_path, 'allmerge.tsv'), sep='\t', index=False)
        unmatched_df.to_csv(os.path.join(sample_path, 'unmatchline.tsv'), sep='\t', index=False)

        print("[INFO] Processed folder: " + sample_dir +
              " | Merged rows: " + str(len(merged_df)) +
              " | Unmatched rows: " + str(len(unmatched_df)))
