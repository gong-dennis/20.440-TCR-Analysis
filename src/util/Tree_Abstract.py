import abc
import pickle as pkl

class Tree_Abstract:
    """
    Abstract class for search trees. Includes methods shared between
    labeled and unlabeled trees, as well as abstract methods for 
    methods that are implemented differently between the child classes.
    """
    def __init__(self, df, root=[], tree_dict={}):
        """
        Constructor for abstract class

        Args:
            df (pd.DataFrame): dataframe containing amino acid data 
                with/without labels
            root (list, optional): list of tree nodes that correspond with the 
                root of the tree (i.e. the first amino acids in the sequence). 
                Defaults to [].
            tree_dict (dict, optional): dictionary storing the tree structure. Datatypes are  Defaults to {}.
        """
        self.df = df
        
        # Can intialize tree without providing data, or with data
        # provided in a dataframe
        if self.df is None:
            self.root = root
            self.tree_dict = tree_dict
        else:
            self.root = []
            self.tree_dict = {}

        # If tree is initialized with data, populate tree from data
        if self.df is not None:
            self._populate_tree()


    def _populate_tree(self):
        """

        """
        for _, row in self.df.iterrows():
            self._add_aa_seq(row)


    @abc.abstractmethod
    def _add_aa_seq(self, row):
        pass


    def save_to_pkl(self, save_file_path):
        with open(save_file_path, 'wb') as f:
            pkl.dump(self, f)


    def save_to_csv(self, save_file_path):
        seqs_df = self._traverse_tree(self)
        seqs_df.to_csv(save_file_path)


    @abc.abstractmethod
    def _traverse_tree(self):
        pass
    
    def get_df(self):
        return self.df
    

    def get_save_file_path(self):
        return self.save_file_path
    

    def get_root(self):
        return self.root
    

    def get_tree_dict(self):
        return self.tree_dict

