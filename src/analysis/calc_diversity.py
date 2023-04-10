import pandas as pd
import os
import numpy as np

def quantify_patients(ov_df, metric_fn, metric_col, ov_new_col):
    # print(f'{metric_col = }')
    new_col_arr = np.zeros(len(ov_df))
    for i, samp in enumerate(ov_df['sample_name']):
        pt_path = os.path.join('data', 'raw', 'ChemoProjTCRs', samp + '.tsv')
        patient_df = pd.read_csv(pt_path, sep='\t', low_memory=False)
        new_col_arr[i] = metric_fn(patient_df, metric_col)
    ov_df[ov_new_col] = new_col_arr
    return ov_df