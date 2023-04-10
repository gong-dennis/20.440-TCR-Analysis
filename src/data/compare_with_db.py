import pandas as pd
import os
from tqdm import tqdm
from src.util.Labeled_AA_Tree import Labeled_AA_Tree
from src.util.Unlabeled_AA_Tree import Unlabeled_AA_Tree
from src.util.helpers import get_aa_seqs

def seq_matches_vdjdb():
    db_path = 'data/raw/vdjdb.txt'
    db_df = pd.read_csv(db_path, sep='\t')
    db_df = db_df[db_df['species'] == 'HomoSapiens']
    db_df = db_df[['cdr3', 'mhc.a', 'mhc.b', 'mhc.class', 'antigen.epitope',	
                  'antigen.gene', 'antigen.species']]
    labeled_tree = Labeled_AA_Tree(db_df)

    save_prefix = 'data/processed/shared_seqs_vdjdb/'
    _generate_matches_db(labeled_tree, save_prefix)


def seq_matches_deepcat_cdr3():
    db_path = 'data/raw/deepcat_cdr3.txt'
    db_df = pd.read_csv(db_path, header=None)
    db_df.rename(columns={0:'cdr3'})
    db_tree = Unlabeled_AA_Tree(db_df)
    save_prefix = 'data/processed/shared_seqs_deepcat/'
    _generate_matches_db(db_tree, save_prefix)


def _generate_matches_db(db_tree, save_prefix):
    data_prefix = 'data/raw/ChemoProjTCRs/'
    file_names = os.listdir(data_prefix)
    for f in tqdm(file_names):
        f_handle = os.path.splitext(f)[0]
        df = _one_patient_matches(db_tree, data_prefix + f)
        df.to_csv(save_prefix + f_handle + '.csv', index=False)


def _one_patient_matches(labeled_tree, livmet_path):
    livmet_df = get_aa_seqs(livmet_path)
    livmet_df = livmet_df.to_frame()
    unlabeled_tree = Unlabeled_AA_Tree(livmet_df)

    lab_comp = labeled_tree.compare_with_unlabeled(unlabeled_tree)
    return lab_comp._traverse_tree()