from src.util.Tree_Abstract import Tree_Abstract
from src.util.Labeled_AA_Tree import Labeled_AA_Tree

import pandas as pd

class Unlabeled_AA_Tree(Tree_Abstract):
    def __init__(self, df, root=[], tree_dict={}):
        super().__init__(df, root, tree_dict)


    def _add_aa_seq(self, row):
        aa_seq = row.iloc[0]
        cum_seq = ''
        for i, aa in enumerate(aa_seq):
            if i == len(aa_seq) - 1:
                current_aa_label = True
            else:
                current_aa_label = False

            cum_seq += aa
            if i == 0:
                if aa not in self.root:
                    self.root.append(aa)
                    self.tree_dict[aa] = [set(), current_aa_label]
            else:
                prev_seq = cum_seq[:-1]
                if cum_seq not in self.tree_dict.keys():
                    self.tree_dict[cum_seq] = [set(), current_aa_label]
                
                if current_aa_label:
                    self.tree_dict[cum_seq][1] = current_aa_label
                self.tree_dict[prev_seq][0].add(cum_seq)


    def compare_with_labeled(self, labeled_tree):
        unlab_root = self.root
        unlab_dict = self.tree_dict
        lab_root = labeled_tree.root
        lab_dict = labeled_tree.tree_dict
        lab_labels = labeled_tree.labels

        shared_roots = list(set(unlab_root).intersection(lab_root))
        shared_nodes = set(unlab_dict.keys()).intersection(lab_dict.keys())
        
        shared_dict = {}
        labels_dict = {}
        for node in shared_nodes:
            unlab_children = unlab_dict[node][0]
            unlab_terminal = unlab_dict[node][1]
            lab_children = lab_dict[node][0]
            lab_terminal = lab_dict[node][1]
            shared_children = unlab_children.intersection(lab_children)
            is_terminal = unlab_terminal and lab_terminal
            shared_dict[node] = (shared_children, is_terminal)
            if is_terminal:
                labels_dict[node] = lab_labels[node]

        return Labeled_AA_Tree(None, root=shared_roots, tree_dict=shared_dict, 
                               labels=labels_dict)


    def compare_with_unlabeled(self, unlabeled_tree):
        unlab1_root = self.root
        unlab1_dict = self.tree_dict
        unlab2_root = unlabeled_tree.root
        unlab2_dict = unlabeled_tree.tree_dict

        shared_roots = list(set(unlab1_root).intersection(unlab2_root))
        shared_nodes = set(unlab1_dict.keys()).intersection(unlab2_dict.keys())
        
        shared_dict = {}
        for node in shared_nodes:
            unlab1_children = unlab1_dict[node][0]
            unlab1_terminal = unlab1_dict[node][1]
            unlab2_children = unlab2_dict[node][0]
            unlab2_terminal = unlab2_dict[node][1]
            shared_children = unlab1_children.intersection(unlab2_children)
            is_terminal = unlab1_terminal and unlab2_terminal
            shared_dict[node] = (shared_children, is_terminal)

        return Unlabeled_AA_Tree(None, root=shared_roots, tree_dict=shared_dict)


    def _traverse_tree(self):
        all_seqs = []
        for aa_seq, val in self.tree_dict.items():
            if val[1]:
                all_seqs.append(aa_seq)
        
        return pd.DataFrame({'aa_seqs': all_seqs}).drop_duplicates()
