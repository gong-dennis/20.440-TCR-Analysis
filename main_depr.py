'''
TODO: 

Make summary CSV files for plotting for:
-- Stripplot for richness and gini impurity --
- Calculate gini from multiplicity of each TCR sequence as well as frequency
appearing in dataset
- Recalculate overview gini values, add overview to summary folder

-- Cluster all datapoints --
- Add all TCRs (without multiplicity) to one array and cluster using clusTCR
- Cross reference clusters with each patient to see which clusters have hits 
for each patient, and see which 



'''



import os
import pandas as pd

import run_parameters as rp

from src.data.calc_overview_metrics import diversity_metrics_to_ov
from src.data.compare_with_db import seq_matches_deepcat_cdr3, seq_matches_vdjdb
from src.data.cluster_tcrs import cluster_patients
from src.data.find_cancer_clusters import find_cancer_clusters
from src.data.tcr_multiplicity import save_tcr_multiplicities
from src.visualization.plot_clustering import survival_vs_deepcat, \
   survival_vs_deepcat_clustering

from src.visualization.plots_from_overview import create_stripplots

def main():
   if rp.compare_with_vdjdb: seq_matches_vdjdb()
   if rp.compare_with_deepcat: seq_matches_deepcat_cdr3()
   if rp.cluster_tcrs: cluster_patients()
   if rp.find_cancer_clusters: find_cancer_clusters()
   if rp.save_tcr_mults: save_tcr_multiplicities()

   saved_df_path = os.path.join('data', 'processed', 'precalculated.tsv')
   if rp.calculate_overview: diversity_metrics_to_ov(saved_df_path)
   ov_df = pd.read_csv(saved_df_path, sep="\t")

   if rp.plot_from_overview: create_stripplots(ov_df)
   if rp.plot_survival_num_tcr: survival_vs_deepcat(ov_df)
   if rp.plot_survival_num_tcr_clustering: survival_vs_deepcat_clustering(ov_df)
    

if __name__ == "__main__":
    main()
