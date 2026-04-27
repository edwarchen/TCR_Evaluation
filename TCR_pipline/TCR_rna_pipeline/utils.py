#!/x03_haplox/users/donglf/miniconda3/bin/python
import os
import sys
import pandas as pd
from collections import Counter
import numpy as np

complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def read_csv(table_file):
    _sep = ',' if table_file.endswith(".csv") else '\t'
    table = pd.read_csv(table_file, sep=_sep)
    return table

def check_file(_file):
    if _file is None or not os.path.isfile(_file):
        print(f"[ERROR] File {_file} does not exist!")
        return False
    return True

def check_duplicates(lst:list, tag:str):
    counts = Counter(lst)
    dup_ids = [k for k,v in counts.items() if v > 1]
    if len(dup_ids) > 0:
        print(f'Duplicates for {tag}')
        print(dup_ids)
    else:
        print(f'No duplication for {tag}')

def check_dir(_dir):
    if not os.path.isdir(_dir):
        print(f"[ERROR] Directory {_dir} does not exist!")
        return False
    return True

def concate_table(dir_lst, table_name):
    _sep = ',' if table_name.endswith(".csv") else '\t'
    table_lst = [pd.read_csv(f"{dir}/{table_name}", sep=_sep) for dir in dir_lst]
    total_table = pd.concat(table_lst, ignore_index=True)
    return total_table

def reverse_seq(seq):
    return seq[::-1].upper()

def reverse_complement(seq):
    return ''.join([complement_dict.get(base, 'N') for base in reverse_seq(seq)])

def gc_content(seq):
    counts = Counter(seq.upper())
    return (counts['G'] + counts['C']) / (counts['G'] + counts['C'] + counts['A'] + counts['T'])


def filter_tbl(input_table:pd.DataFrame, count_threshold, len_low, len_up) -> pd.DataFrame:
    """
    filter input tcr seq tbl
    """
    input_table = input_table[input_table['count'] >= count_threshold].copy() # avoid warning information
    input_table['cdr3_len'] = input_table['cdr3aa'].map(len)
    input_table = input_table[(input_table['cdr3_len'] <= len_up) & (input_table['cdr3_len'] >= len_low)]
    input_table = input_table[~input_table['cdr3aa'].str.contains('[_\*]')] # filter cdr3aa that got '_' or '*' in it.
    input_table = input_table[input_table['cdr3aa'].str.startswith('C')]
    input_table = input_table[input_table['cdr3aa'].str.endswith('F')]
    total_count = input_table['count'].sum()
    input_table['cdr3aa_freq'] = input_table['count'] / total_count
    input_table['score'] = input_table['count'].map(np.log)
    return input_table


def filter_tbl_nfunc(input_table:pd.DataFrame) -> pd.DataFrame:
    """
    filter input tcr seq tbl
    """
    input_table = input_table[~input_table['cdr3aa'].str.contains('[_\*]')] # filter cdr3aa that got '_' or '*' in it.
    input_table = input_table[input_table['cdr3aa'].str.startswith('C')]
    input_table = input_table[input_table['cdr3aa'].str.endswith('F')]
    input_table['cdr3_len'] = input_table['cdr3aa'].map(len)
    total_count = input_table['count'].sum()
    input_table['cdr3aa_freq'] = input_table['count'] / total_count
    input_table['score'] = input_table['count'].map(np.log)
    return input_table

def filter_tbl_threshold(input_table:pd.DataFrame, count_threshold:int) -> pd.DataFrame:
    """
    filter input tcr seq tbl
    """
    input_table = input_table[input_table['count'] >= count_threshold].copy() # avoid warning information
    return input_table

def commpare_tbl(tab1, tab2, id_col):
    col1 = set(tab1.columns.tolist())
    col2 = set(tab2.columns.tolist())
    common_cols = col1 & col2
    sort_tab1 = tab1.sort_values(by=[id_col])
    sort_tab2 = tab2.sort_values(by=[id_col])
    print(common_cols)
    print(col1 - common_cols)
    print(col2 - common_cols)
    for _col in common_cols:
        if sort_tab1[_col].to_list() != sort_tab2[_col].tolist():
            print(_col)