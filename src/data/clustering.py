from os import listdir
from os.path import splitext, join
import pandas as pd
from clustcr import Clustering
from tqdm import tqdm

def cluster_each_patient_seq(seq_condensed_dir, each_pt_cluster_dir):
    file_names = listdir(seq_condensed_dir)
    num_pts = len(file_names)
    
    for i, f in enumerate(file_names):
        print(f' - Progress: {i+1} of {num_pts}')
        df = pd.read_csv(join(seq_condensed_dir, f), sep='\t', low_memory=False)
        aa_series = df['aminoAcid']
        clst = Clustering(method='mcl')
        output = clst.fit(aa_series)
        output.clusters_df.to_csv(join(each_pt_cluster_dir, f), index=False)


def cluster_all_patient_seqs(seq_condensed_dir, all_pt_cluster_path):
    file_names = listdir(seq_condensed_dir)
    num_pts = len(file_names)
    
    all_aas = set()
    for i, f in enumerate(file_names):
        print(f' - Progress: {i+1} of {num_pts}')
        df = pd.read_csv(join(seq_condensed_dir, f), sep='\t', low_memory=False)
        all_aas.update(df['aminoAcid'])
        
    clst = Clustering(method='two-step')
    output = clst.fit(pd.Series(list(all_aas)))
    output.clusters_df.to_csv(all_pt_cluster_path, index=False)


def add_cancer_assoc_annot(deepcat_db_path, all_pt_cluster_path):
    superc_df = pd.read_csv(all_pt_cluster_path, sep=',', low_memory=False)
    db_df = db_df = pd.read_csv(deepcat_db_path, header=None)
    db_set = set(db_df[0])

    cancer_clust_set = set()
    for i, row in tqdm(superc_df.iterrows()):
        if row['junction_aa'] in db_set:
            cancer_clust_set.add(row['cluster'])
    
    is_assoc = [c in cancer_clust_set for c in superc_df['cluster']]
    superc_df['is_cancer_assoc'] = is_assoc
    superc_df.to_csv(all_pt_cluster_path, index=False, sep=',')