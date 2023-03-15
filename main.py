import os

from src.data.parse_overview import parse_sample_overview
from src.analysis.calc_diversity import quantify_patients
from src.util.helpers import get_mean, shannon_idx, shannon_idx_vdj, \
       get_richness, get_richness_vdj
from src.visualization.plotting import violin

def main():
    load_precalculated = False
    
    if not load_precalculated:
       # Parse overview tsv file
       samp_overview_path = './data/analysis/SampleOverview.tsv'
       ov_df = parse_sample_overview(samp_overview_path)

       # Calculate from patient sequencing data
       metric_fns_cols = [(shannon_idx, 'aminoAcid', 'aa_shannon'),
                            (get_mean, 'cdr3Length', 'mean_cdr3'),
                            (shannon_idx, 'vGeneName', 'shannon_vgene'),
                            (shannon_idx, 'jGeneName', 'shannon_jgene'),
                            (shannon_idx_vdj, 'vdjGenes', 'shannon_vdjgene'),
                            (get_richness, 'aminoAcid', 'aa_richness'),
                            (get_richness, 'vGeneName', 'richness_vgene'),
                            (get_richness, 'jGeneName', 'richness_jgene'),
                            (get_richness_vdj, 'vdjGenes', 'richness_vdjgene'),]
       for metric_tup in metric_fns_cols:
              ov_df = quantify_patients(ov_df, 
                                   metric_tup[0], metric_tup[1], metric_tup[2])
       ov_df.to_csv(os.path.join('data', 'processed', 'precalculated.tsv'), 
                    sep="\t")
    
    ## Plotting
    xlabel = 'Patient Group'
    violin(ov_df, 'interval_group_num', 'fraction_productive', 
           xlabel, 'Fraction Productive', 'fraction_productive_violin.jpg')
    violin(ov_df, 'interval_group_num', 'rel_prod_rearr', 
           xlabel, 'Relative Productive Rearrangements', 'rel_prod_rearr_violin.jpg')
    violin(ov_df, 'interval_group_num', 'rel_total_rearr', 
           xlabel, 'Relative Total Rearrangements', 'rel_total_rearr_violin.jpg')
    violin(ov_df, 'interval_group_num', 'frac_prod_rearr', 
           xlabel, 'Fraction of Productive Rearrangements', 'frac_prod_rearr_violin.jpg')
    violin(ov_df, 'interval_group_num', 'productive_simpson_clonality', 
           xlabel, 'Productive Simpson Clonality', 'prod_simps_clon_rearr_violin.jpg')
    violin(ov_df, 'interval_group_num', 'max_productive_frequency', 
           xlabel, 'Max Productive Frequency', 'max_prod_freq_violin.jpg')
    
    violin(ov_df, 'interval_group_num', 'aa_shannon', 
           xlabel, 'Shannon Index (from AA seqs)', 'aa_shannon_violin.jpg')
    violin(ov_df, 'interval_group_num', 'mean_cdr3', 
           xlabel, 'Mean CDR3 Length', 'mean_cdr3_violin.jpg')
    violin(ov_df, 'interval_group_num', 'shannon_vgene', 
           xlabel, 'Shannon Index (TCR V-gene)', 'shannon_vgene_violin.jpg')
    violin(ov_df, 'interval_group_num', 'shannon_jgene', 
           xlabel, 'Shannon Index (TCR J-gene)', 'shannon_jgene_violin.jpg')
    violin(ov_df, 'interval_group_num', 'shannon_vdjgene', 
           xlabel, 'Shannon Index (from V, D, J seqs)', 'shannon_vdjgene_violin.jpg')
    
    violin(ov_df, 'interval_group_num', 'aa_richness', 
           xlabel, 'Richness (from AA seqs)', 'rich_shannon_violin.jpg')
    violin(ov_df, 'interval_group_num', 'richness_vgene', 
           xlabel, 'Richness (TCR V-gene)', 'rich_vgene_violin.jpg')
    violin(ov_df, 'interval_group_num', 'richness_jgene', 
           xlabel, 'Richness (TCR J-gene)', 'rich_jgene_violin.jpg')
    violin(ov_df, 'interval_group_num', 'richness_vdjgene', 
           xlabel, 'Richness (from V, D, J seqs)', 'rich_vdjgene_violin.jpg')
    


if __name__ == "__main__":
    main()