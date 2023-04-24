from src.analysis.calc_diversity import quantify_patients
from src.util.pop_metric_helpers import get_richness_vdj, gini_impurity_vdj, \
    shannon_idx_vdj, gini_tcr_mults, shannon_tcr_mults

from os import listdir
from os.path import join, splitext
from tqdm import tqdm
import pandas as pd
import numpy as np
from scipy.stats import rankdata

def condense_seq_data(raw_patient_dir, seq_condensed_dir):
    file_names = listdir(raw_patient_dir)
    for f in tqdm(file_names):
        full_path = join(raw_patient_dir, f)
        f_handle = splitext(f)[0]
        df = pd.read_csv(full_path, sep="\t", low_memory=False)
        
        aa_codes = 'ARNDCEQGHILKMFPSTWYVX'
        df = df.loc[df['aminoAcid'].notna()]
        df = df.loc[df['aminoAcid'].apply(
            lambda x: all(i in aa_codes for i in x))]

        output_dict = {c: [] for c in df.columns.to_list()}
        for _, row in df.iterrows():
            if row['aminoAcid'] not in output_dict['aminoAcid']:
                for k in output_dict.keys():
                    output_dict[k].append(row[k])
            else:
                list_idx = output_dict['aminoAcid'].index(row['aminoAcid'])
                output_dict['count (templates/reads)'][list_idx] += \
                    row['count (templates/reads)']
        
        output = pd.DataFrame(output_dict)
        output.to_csv(join(seq_condensed_dir, f_handle + '.csv'), index=False, 
                      sep="\t")

def add_overview_metrics(ov_parsed_path, seq_condensed_dir, added_metrics_path):
    file_names = listdir(seq_condensed_dir)
    ordered_idx = rankdata([int(i[7:9]) for i in file_names]).astype(int) - 1
    
    richness_aa = np.zeros(len(file_names))
    richness_vdj = np.zeros(len(file_names))
    gini_aa = np.zeros(len(file_names))
    gini_vdj = np.zeros(len(file_names))
    shannon_aa = np.zeros(len(file_names))
    shannon_vdj = np.zeros(len(file_names))
    
    for i, f in enumerate(tqdm(file_names)):
        df = pd.read_csv(join(seq_condensed_dir, f), sep="\t", low_memory=False)
        corr_i = ordered_idx[i]

        richness_aa[corr_i] = len(df)
        richness_vdj[corr_i] = get_richness_vdj(df, 'vdjGenes')
        
        gini_aa[corr_i] = gini_tcr_mults(df)
        gini_vdj[corr_i] = gini_impurity_vdj(df, 'vdjGenes')

        shannon_aa[corr_i] = shannon_tcr_mults(df)
        shannon_vdj[corr_i] = shannon_idx_vdj(df, 'vdjGenes')
    
    output = pd.read_csv(ov_parsed_path, sep="\t", low_memory=False)
    output['richness_aa'] = richness_aa
    output['richness_vdjgene'] = richness_vdj
    output['gini_aa'] = gini_aa
    output['gini_vdjgene'] = gini_vdj
    output['shannon_aa'] = shannon_aa
    output['shannon_vdjgene'] = shannon_vdj
    output.to_csv(added_metrics_path, index=False, 
                  sep="\t")