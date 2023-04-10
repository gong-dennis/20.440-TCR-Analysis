from clustcr import Clustering

import os
import pandas as pd
from src.util.helpers import get_aa_seqs
import matplotlib.pyplot as plt

def cluster_patients():
    save_prefix = 'data/processed/clustcr_labels/'
    data_prefix = 'data/raw/ChemoProjTCRs/'
    file_names = os.listdir(data_prefix)
    
    print('Clustering patient data...')
    for f in file_names:
        livmet_df = get_aa_seqs(data_prefix + f).drop_duplicates()
        f_handle = os.path.splitext(f)[0]
        clst = Clustering(method='mcl')
        output = clst.fit(livmet_df)
        output.clusters_df.to_csv(save_prefix + f_handle + '.csv', index=False)