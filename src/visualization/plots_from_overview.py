from src.visualization.plotting import strip_plot

def create_stripplots(ov_df):
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