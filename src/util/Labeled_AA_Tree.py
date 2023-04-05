from src.util.Tree_Abstract import Tree_Abstract

import pandas as pd

class Labeled_AA_Tree(Tree_Abstract): 
    def __init__(self, df, root=[], tree_dict={}, labels={}):
        self.labels = labels
        super().__init__(df, root, tree_dict)

    
    def _add_aa_seq(self, row):
        aa_seq = row.iloc[0]
        label = row.iloc[1:]
        cum_seq = ''
        for i, aa in enumerate(aa_seq):
            cum_seq += aa
            if i == len(aa_seq) - 1:
                current_aa_label = True
                if cum_seq in self.labels.keys():
                    self.labels[cum_seq].append(label)
                else:
                    self.labels[cum_seq] = [label]
            else:
                current_aa_label = False

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


    def compare_with_unlabeled(self, unlabeled_tree):
        unlab_root = unlabeled_tree.root
        unlab_dict = unlabeled_tree.tree_dict
        lab_root = self.root
        lab_dict = self.tree_dict
        lab_labels = self.labels

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
    

    def _traverse_tree(self):
        all_seqs = []
        seq_labels = []
        for aa_seq, val in self.tree_dict.items():
            is_terminal = val[1]
            if is_terminal:
                for lbl in self.labels[aa_seq]:
                    all_seqs.append(aa_seq)
                    seq_labels.append(lbl) # THIS IS TERRIBLE

        return pd.DataFrame({'aa_seqs': all_seqs, 'labels': seq_labels})