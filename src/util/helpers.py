import numpy as np
from collections import Counter


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


def get_richness(patient_df, col_name):
    return len(patient_df[col_name].unique())


def get_richness_vdj(patient_df, col_name):
    patient_df[col_name] = patient_df['vGeneName'] + ', ' \
        + patient_df['dGeneName'] + ', ' + patient_df['jGeneName']
    return get_richness(patient_df, col_name)