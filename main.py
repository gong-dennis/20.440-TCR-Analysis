import os
import pandas as pd

import run_parameters as rp

from src.data.parse_overview import parse_sample_overview
from src.analysis.calc_diversity import quantify_patients
from src.util.helpers import get_mean, shannon_idx, shannon_idx_vdj, \
       get_richness, get_richness_vdj, gini_impurity, gini_impurity_vdj
from src.data.compare_with_db import seq_matches_deepcat_cdr3, seq_matches_vdjdb

from src.visualization.plots_from_overview import create_stripplots

def main():
<<<<<<< Updated upstream
    if rp.compare_with_vdjdb: seq_matches_vdjdb()
    if rp.compare_with_deepcat: seq_matches_deepcat_cdr3()

=======
    load_precalculated = False
    
>>>>>>> Stashed changes
    saved_df_path = os.path.join('data', 'processed', 'precalculated.tsv')
    if not rp.load_overview:
       # Parse overview tsv file
       samp_overview_path = './data/analysis/SampleOverview.tsv'
       ov_df = parse_sample_overview(samp_overview_path)

       # Calculate from patient sequencing data
       metric_fns_cols = [(shannon_idx, 'aminoAcid', 'aa_shannon'),
                          (get_mean, 'cdr3Length', 'mean_cdr3'),
                          (shannon_idx, 'vGeneName', 'shannon_vgene'),
                          (shannon_idx, 'jGeneName', 'shannon_jgene'),
                          (shannon_idx_vdj, 'vdjGenes', 'shannon_vdjgene'),
                          (gini_impurity, 'vGeneName', 'shannon_vgene'),
                          (gini_impurity, 'jGeneName', 'shannon_jgene'),
                          (gini_impurity_vdj, 'vdjGenes', 'shannon_vdjgene'),
                          (get_richness, 'aminoAcid', 'aa_richness'),
                          (get_richness, 'vGeneName', 'richness_vgene'),
                          (get_richness, 'jGeneName', 'richness_jgene'),
                          (get_richness_vdj, 'vdjGenes', 'richness_vdjgene')]
       for metric_tup in metric_fns_cols:
              ov_df = quantify_patients(ov_df, 
                                   metric_tup[0], metric_tup[1], metric_tup[2])
       ov_df.to_csv(saved_df_path, sep="\t")
    else:
       ov_df = pd.read_csv(saved_df_path, sep="\t")

    if rp.plot_from_overview: create_stripplots(ov_df)
    

if __name__ == "__main__":
    main()
