import pandas as pd
import numpy as np
import os

def generate_shared_seqs_db():
    db_path = 'data/raw/vdjdb.txt'
    db_df = pd.read_csv(db_path, sep='\t')
    db_df = db_df[db_df['species'] == 'HomoSapiens']
    db_df = db_df[['cdr3', 'mhc.a', 'mhc.b', 'mhc.class', 'antigen.epitope',	
                'antigen.gene', 'antigen.species']]
    
    data_prefix = 'data/raw/ChemoProjTCRs/'
    save_prefix = 'data/processed/shared_seqs_vdjdb'
    file_names = os.listdir('data/raw/ChemoProjTCRs')
    # file_names = ['LivMet_35.tsv', 'LivMet_89.tsv', 'LivMet_01.tsv']
    for i, f in enumerate(file_names):
        print(f'iter = {i+1} of {len(file_names)}')
        f_handle = os.path.splitext(f)[0]
        df = _generate_one_sequence(db_df, data_prefix + f)
        df.to_csv(save_prefix + f_handle + '.csv', index=False)


def _generate_one_sequence(db_df, livmet_path):
    aa_codes = 'ARNDCEQGHILKMFPSTWYVX'

    livmet_df = pd.read_csv(livmet_path, sep='\t')
    livmet_df = livmet_df.loc[livmet_df['aminoAcid'].notna()]
    livmet_df = livmet_df.loc[livmet_df['aminoAcid'].apply(
        lambda x: all(i in aa_codes for i in x))]
    livmet_df = livmet_df['aminoAcid'].to_frame()

    unlab_list = livmet_df['aminoAcid'].to_list()
    lab_list = db_df['cdr3'].to_list()
    sl = {}
    for i in unlab_list:
        if i in lab_list:
            slice = db_df.loc[db_df['cdr3']==i]
            keys = slice.columns.to_list()
            if not sl:
                sl = {i: [] for i in keys}

            for _, row in slice.iterrows():
                for j in keys:
                    sl[j].append(row[j])
    return pd.DataFrame(sl).drop_duplicates()