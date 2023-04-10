import os
from tqdm import tqdm
import pandas as pd
from src.util.Labeled_AA_Tree import Labeled_AA_Tree
from src.util.Unlabeled_AA_Tree import Unlabeled_AA_Tree

def find_cancer_clusters():
    save_prefix = 'data/processed/cancer_clusters/'
    data_prefix = 'data/processed/clustcr_labels/'
    file_names = os.listdir(data_prefix)
    for f in tqdm(file_names):
        cancer_assoc = pd.read_csv('data/processed/shared_seqs_deepcat/' + f)
        cancer_tree = Unlabeled_AA_Tree(cancer_assoc)

        clusters = pd.read_csv(data_prefix + f)
        cluster_tree = Labeled_AA_Tree(clusters)

        comp = cluster_tree.compare_with_unlabeled(cancer_tree)
        comp._traverse_tree().to_csv(save_prefix + f, index=False)