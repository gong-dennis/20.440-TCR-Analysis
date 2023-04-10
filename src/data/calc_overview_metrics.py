from tqdm import tqdm

from src.data.parse_overview import parse_sample_overview
from src.analysis.calc_diversity import quantify_patients
from src.util.helpers import get_mean, shannon_idx, shannon_idx_vdj, \
       get_richness, get_richness_vdj, gini_impurity, gini_impurity_vdj

def diversity_metrics_to_ov(saved_df_path):
    # Parse overview tsv file
       samp_overview_path = './data/analysis/SampleOverview.tsv'
       ov_df = parse_sample_overview(samp_overview_path)

       # Calculate from patient sequencing data
       metric_fns_cols = [(shannon_idx, 'aminoAcid', 'aa_shannon'),
                          (get_mean, 'cdr3Length', 'mean_cdr3'),
                          (shannon_idx, 'vGeneName', 'shannon_vgene'),
                          (shannon_idx, 'jGeneName', 'shannon_jgene'),
                          (shannon_idx_vdj, 'vdjGenes', 'shannon_vdjgene'),
                          (gini_impurity, 'aminoAcid', 'aa_gini'),
                          (gini_impurity, 'vGeneName', 'gini_vgene'),
                          (gini_impurity, 'jGeneName', 'gini_jgene'),
                          (gini_impurity_vdj, 'vdjGenes', 'gini_vdjgene'),
                          (get_richness, 'aminoAcid', 'aa_richness'),
                          (get_richness, 'vGeneName', 'richness_vgene'),
                          (get_richness, 'jGeneName', 'richness_jgene'),
                          (get_richness_vdj, 'vdjGenes', 'richness_vdjgene')]
       
       print('Calculating additional metrics for overview file...')
       for metric_tup in tqdm(metric_fns_cols):
              ov_df = quantify_patients(ov_df, 
                                   metric_tup[0], metric_tup[1], metric_tup[2])
       ov_df.to_csv(saved_df_path, sep="\t")