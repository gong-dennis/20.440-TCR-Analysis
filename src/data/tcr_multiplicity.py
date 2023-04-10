from src.util.helpers import get_aa_seqs
import os
from tqdm import tqdm
import pickle as pkl

def save_tcr_multiplicities():
    data_prefix = 'data/raw/ChemoProjTCRs/'
    file_names = os.listdir(data_prefix)
    for f in tqdm(file_names):
        f_handle = os.path.splitext(f)[0]
        aa_df = get_aa_seqs(data_prefix + f, whole_df=True)
        mult_dict = {}
        for _, row in aa_df.iterrows():
            if row['aminoAcid'] in mult_dict.keys():
                mult_dict[row['aminoAcid']] += row['count (templates/reads)']
            else:
                mult_dict[row['aminoAcid']] = row['count (templates/reads)']
        pkl_path = 'data/processed/multiplicity_dicts/'
        pkl.dump(mult_dict, open(pkl_path + f_handle + '.p', 'wb'))


def get_tcr_multiplicity(path):
    return pkl.load(open(path, 'rb'))