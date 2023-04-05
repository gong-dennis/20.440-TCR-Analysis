from src.util.Labeled_AA_Tree import Labeled_AA_Tree
from src.util.Unlabeled_AA_Tree import Unlabeled_AA_Tree

import pandas as pd
import numpy as np

# labeled_seqs = ['abcd', 'acd', 'abb', 'abcde', 'acb', 'ac']
# labeled_labels = ['Human', 'HIV', 'Human', 'HPV', 'HIV', 'HIV']
# labeled_seq_df = pd.DataFrame({'aa_df': labeled_seqs, 'labels': labeled_labels})

# unlabeled_seqs = ['acb', 'abcd', 'ab']
# unlabeled_seq_df = pd.DataFrame({'aa_df': unlabeled_seqs})

# unlabeled_seqs2 = ['acb', 'abcde', 'ab', 'a']
# unlabeled_seq2_df = pd.DataFrame({'aa_df': unlabeled_seqs2})

# labeled_tree = Labeled_AA_Tree(labeled_seq_df)
# unlabeled_tree = Unlabeled_AA_Tree(unlabeled_seq_df)
# unlabeled_tree2 = Unlabeled_AA_Tree(unlabeled_seq2_df)

# l1 = labeled_tree.compare_with_unlabeled(unlabeled_tree)
# print(l1._traverse_tree())
# l2 = unlabeled_tree.compare_with_labeled(labeled_tree)
# print(l2._traverse_tree())
# l3 = unlabeled_tree.compare_with_unlabeled(unlabeled_tree2)
# print(l3._traverse_tree())

aa_codes = 'ARNDCEQGHILKMFPSTWYVX'

livmet_path = 'data/raw/ChemoProjTCRs/LivMet_01.tsv'
livmet_df = pd.read_csv(livmet_path, sep='\t')
livmet_df = livmet_df.loc[livmet_df['aminoAcid'].notna()]
livmet_df = livmet_df.loc[livmet_df['aminoAcid'].apply(lambda x: all(i in aa_codes for i in x))]
livmet_df = livmet_df['aminoAcid'].to_frame()
unlabeled_tree = Unlabeled_AA_Tree(livmet_df)

db_path = 'data/raw/vdjdb.txt'
db_df = pd.read_csv(db_path, sep='\t')
db_df = db_df[db_df['species'] == 'HomoSapiens']
db_df = db_df[['cdr3', 'mhc.a', 'mhc.b', 'mhc.class', 'antigen.epitope',	
              'antigen.gene', 'antigen.species']]
labeled_tree = Labeled_AA_Tree(db_df)
'a'


# df['aminoAcid'].apply(lambda x: 'ASS' in x)
# def filter_aa_chars(series):
    