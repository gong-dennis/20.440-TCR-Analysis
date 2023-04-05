import os
import pandas as pd

from src.data.parse_overview import parse_sample_overview
from src.analysis.calc_diversity import quantify_patients
from src.util.helpers import get_mean, shannon_idx, shannon_idx_vdj, \
       get_richness, get_richness_vdj
from src.visualization.plotting import strip_plot

def main():
    load_precalculated = True
    
    saved_df_path = os.path.join('data', 'processed', 'precalculated.tsv')
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
                          (get_richness_vdj, 'vdjGenes', 'richness_vdjgene')]
       for metric_tup in metric_fns_cols:
              ov_df = quantify_patients(ov_df, 
                                   metric_tup[0], metric_tup[1], metric_tup[2])
       ov_df.to_csv(saved_df_path, sep="\t")
    else:
       ov_df = pd.read_csv(saved_df_path, sep="\t")

    ## Plotting
    xslice = 'group_label'
    strip_plot(ov_df, xslice, 'fraction_productive', 
           'Fraction Productive', 'fraction_productive_violin.jpg')
    strip_plot(ov_df, xslice, 'rel_prod_rearr', 
           'Relative Productive Rearrangements', 'rel_prod_rearr_violin.jpg')
    strip_plot(ov_df, xslice, 'rel_total_rearr', 
           'Relative Total Rearrangements', 'rel_total_rearr_violin.jpg')
    strip_plot(ov_df, xslice, 'frac_prod_rearr', 
           'Fraction of Productive Rearrangements', 'frac_prod_rearr_violin.jpg')
    strip_plot(ov_df, xslice, 'productive_simpson_clonality', 
           'Productive Simpson Clonality', 'prod_simps_clon_rearr_violin.jpg')
    strip_plot(ov_df, xslice, 'max_productive_frequency', 
           'Max Productive Frequency', 'max_prod_freq_violin.jpg')
    
    strip_plot(ov_df, xslice, 'aa_shannon', 
           'Shannon Index (from Amino Acid Sequences)', 'aa_shannon_violin.jpg')
    strip_plot(ov_df, xslice, 'mean_cdr3', 
           'Mean CDR3 Length', 'mean_cdr3_violin.jpg')
    strip_plot(ov_df, xslice, 'shannon_vgene', 
           'Shannon Index (TCR V-gene)', 'shannon_vgene_violin.jpg')
    strip_plot(ov_df, xslice, 'shannon_jgene', 
           'Shannon Index (TCR J-gene)', 'shannon_jgene_violin.jpg')
    strip_plot(ov_df, xslice, 'shannon_vdjgene', 
           'Shannon Index (Unique V, D, J-gene Combinations)', 
           'shannon_vdjgene_violin.jpg')
    
    strip_plot(ov_df, xslice, 'aa_richness', 
           'Richness (from Amino Acid Sequences)', 'rich_shannon_violin.jpg')
    strip_plot(ov_df, xslice, 'richness_vgene', 
           'Richness (TCR V-gene)', 'rich_vgene_violin.jpg')
    strip_plot(ov_df, xslice, 'richness_jgene', 
           'Richness (TCR J-gene)', 'rich_jgene_violin.jpg')
    strip_plot(ov_df, xslice, 'richness_vdjgene', 
           'Richness (Unique V, D, J-gene Combinations)', 
           'rich_vdjgene_violin.jpg')
    

if __name__ == "__main__":
    main()