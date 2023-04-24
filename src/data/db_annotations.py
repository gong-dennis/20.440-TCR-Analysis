from os import listdir
from os.path import join
import pandas as pd
from tqdm import tqdm

def get_deepcat_annotations(deepcat_db_path, seq_condensed_dir, 
                            deepcat_annotate_dir):
    db_df = pd.read_csv(deepcat_db_path, header=None)
    db_df = db_df.rename(columns={0:'cdr3'})
    db_set = set(db_df['cdr3'])

    file_names = listdir(seq_condensed_dir)
    
    for _, f in enumerate(tqdm(file_names)):
        df = pd.read_csv(join(seq_condensed_dir, f), sep='\t', low_memory=False)
        patient_set = set(df['aminoAcid'])
        shared_seqs = db_set.intersection(patient_set)
        save_df = pd.DataFrame({'aminoAcid': list(shared_seqs)})
        save_df.to_csv(join(deepcat_annotate_dir, f), index=False)


def annotate_overview(ov_added_metrics_path, deepcat_annotate_dir, 
                         each_pt_cluster_dir, seq_condensed_dir, 
                         ov_num_cancer_tcrs_path):
    pt_seqs_names = sorted(listdir(seq_condensed_dir))

    cum_list = []
    for _, f in enumerate(tqdm(pt_seqs_names)):
        cluster_file_path = join(each_pt_cluster_dir, f)
        annotation_file_path = join(deepcat_annotate_dir, f)
        pt_seq_file_path = join(seq_condensed_dir, f)

        clust_df = pd.read_csv(cluster_file_path, sep=',', low_memory=False)
        anno_df = pd.read_csv(annotation_file_path, sep='\t', low_memory=False)
        seq_df = pd.read_csv(pt_seq_file_path, sep='\t', low_memory=False)

        cancer_tcrs = anno_df.iloc[:,0].to_list()
        cum_mult = 0
        for seq in cancer_tcrs:
            if seq in clust_df['junction_aa'].values:
                slice = clust_df.loc[clust_df['junction_aa'] == seq]
                cluster_num = slice['cluster'].iloc[0]
                cluster_members = clust_df.loc[
                    clust_df['cluster'] == cluster_num]['junction_aa'].to_list()
                for c in cluster_members:
                    # cum_mult += mults_dict[c]
                    cum_mult += _find_mult(seq_df, c)
            else:
                # cum_mult += mults_dict[seq]
                cum_mult += _find_mult(seq_df, seq)
        cum_list.append(cum_mult)
    
    ov_df = pd.read_csv(ov_added_metrics_path, sep='\t', low_memory=False)
    ov_df['num_cancer_tcrs'] = cum_list
    ov_df.to_csv(ov_num_cancer_tcrs_path, index=False, sep='\t')

def _find_mult(seq_df, c):
    num_mults = seq_df.loc[seq_df['aminoAcid'] == c, 
                           'count (templates/reads)'].values[0]
    return num_mults