import numpy as np
from collections import Counter
import pandas as pd

def get_mean(patient_df, col_name):
    return np.mean(patient_df[col_name])


def shannon_idx(patient_df, col_name):
    counts = Counter(patient_df[col_name])
    if np.NaN in counts.keys():
        del counts[np.NaN]
    total_count = sum(counts.values())
    freqs = np.array([v / total_count for _, v in counts.items()])
    return -np.sum(freqs * np.log(freqs))


def shannon_idx_vdj(patient_df, col_name):
    patient_df[col_name] = patient_df['vGeneName'] + ', ' \
        + patient_df['dGeneName'] + ', ' + patient_df['jGeneName']
    return shannon_idx(patient_df, col_name)


def gini_impurity(patient_df, col_name):
    counts = Counter(patient_df[col_name])
    if np.NaN in counts.keys():
        del counts[np.NaN]
    total_count = sum(counts.values())
    freqs = np.array([v / total_count for _, v in counts.items()])
    return 1 - np.sum(freqs ** 2)


def gini_impurity_vdj(patient_df, col_name):
    patient_df[col_name] = patient_df['vGeneName'] + ', ' \
        + patient_df['dGeneName'] + ', ' + patient_df['jGeneName']
    return gini_impurity(patient_df, col_name)


def get_richness(patient_df, col_name):
    return len(patient_df[col_name].unique())


def get_richness_vdj(patient_df, col_name):
    patient_df[col_name] = patient_df['vGeneName'] + ', ' \
        + patient_df['dGeneName'] + ', ' + patient_df['jGeneName']
    return get_richness(patient_df, col_name)


def get_aa_seqs(livmet_path, whole_df=False):
    aa_codes = 'ARNDCEQGHILKMFPSTWYVX'

    livmet_df = pd.read_csv(livmet_path, sep='\t', low_memory=False)
    livmet_df = livmet_df.loc[livmet_df['aminoAcid'].notna()]
    livmet_df = livmet_df.loc[livmet_df['aminoAcid'].apply(
        lambda x: all(i in aa_codes for i in x))]
    
    if whole_df:
        return livmet_df
    else:
        return livmet_df['aminoAcid']