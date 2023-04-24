import numpy as np
from collections import Counter
import pandas as pd

def gini_tcr_mults(patient_df):
    total = np.sum(patient_df['count (templates/reads)'])
    freqs = np.array([v / total for v in patient_df['count (templates/reads)']])
    return 1 - np.sum(freqs ** 2)

def shannon_tcr_mults(patient_df):
    total = np.sum(patient_df['count (templates/reads)'])
    freqs = np.array([v / total for v in patient_df['count (templates/reads)']])
    return -np.sum(freqs * np.log(freqs))

def gini_unique_vals(patient_df, col_name):
    counts = Counter(patient_df[col_name])
    if np.NaN in counts.keys():
        del counts[np.NaN]
    total_count = sum(counts.values())

    patient_df[col_name].to_list()
    freqs = np.array([v / total_count for _, v in counts.items()])
    return 1 - np.sum(freqs ** 2)


def gini_impurity_vdj(patient_df, col_name):
    patient_df[col_name] = patient_df['vGeneName'] + ', ' \
        + patient_df['dGeneName'] + ', ' + patient_df['jGeneName']
    return gini_unique_vals(patient_df, col_name)


def shannon_unique_vals(patient_df, col_name):
    counts = Counter(patient_df[col_name])
    if np.NaN in counts.keys():
        del counts[np.NaN]
    total_count = sum(counts.values())
    freqs = np.array([v / total_count for _, v in counts.items()])
    return -np.sum(freqs * np.log(freqs))


def shannon_idx_vdj(patient_df, col_name):
    patient_df[col_name] = patient_df['vGeneName'] + ', ' \
        + patient_df['dGeneName'] + ', ' + patient_df['jGeneName']
    return shannon_unique_vals(patient_df, col_name)


def get_richness(patient_df, col_name):
    return len(patient_df[col_name].unique())


def get_richness_vdj(patient_df, col_name):
    patient_df[col_name] = patient_df['vGeneName'] + ', ' \
        + patient_df['dGeneName'] + ', ' + patient_df['jGeneName']
    return get_richness(patient_df, col_name)