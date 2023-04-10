from src.visualization.plotting import strip_plot

def create_stripplots(ov_df):
    ## Plotting
    xslice = 'group_label'
    strip_plot(ov_df, xslice, 'fraction_productive', 
           'Fraction Productive', 'fraction_productive_strip.jpg')
    strip_plot(ov_df, xslice, 'rel_prod_rearr', 
           'Relative Productive Rearrangements', 'rel_prod_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'rel_total_rearr', 
           'Relative Total Rearrangements', 'rel_total_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'frac_prod_rearr', 
           'Fraction of Productive Rearrangements', 'frac_prod_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'productive_simpson_clonality', 
           'Productive Simpson Clonality', 'prod_simps_clon_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'max_productive_frequency', 
           'Max Productive Frequency', 'max_prod_freq_strip.jpg')
    
    strip_plot(ov_df, xslice, 'aa_shannon', 
           'Shannon Index (from Amino Acid Sequences)', 'aa_shannon_strip.jpg')
    strip_plot(ov_df, xslice, 'mean_cdr3', 
           'Mean CDR3 Length', 'mean_cdr3_strip.jpg')
    strip_plot(ov_df, xslice, 'shannon_vgene', 
           'Shannon Index (TCR V-gene)', 'shannon_vgene_strip.jpg')
    strip_plot(ov_df, xslice, 'shannon_jgene', 
           'Shannon Index (TCR J-gene)', 'shannon_jgene_strip.jpg')
    strip_plot(ov_df, xslice, 'shannon_vdjgene', 
           'Shannon Index (Unique V, D, J-gene Combinations)', 
           'shannon_vdjgene_strip.jpg')
    
    strip_plot(ov_df, xslice, 'aa_gini', 
           'Gini Impurity (from Amino Acid Sequences)', 'aa_gini_strip.jpg')
    strip_plot(ov_df, xslice, 'gini_vgene', 
           'Gini Impurity (TCR V-gene)', 'gini_vgene_strip.jpg')
    strip_plot(ov_df, xslice, 'gini_jgene', 
           'Gini Impurity (TCR J-gene)', 'gini_jgene_strip.jpg')
    strip_plot(ov_df, xslice, 'gini_vdjgene', 
           'Gini Impurity (Unique V, D, J-gene Combinations)', 
           'gini_vdjgene_strip.jpg')
    
    strip_plot(ov_df, xslice, 'aa_richness', 
           'Richness (from Amino Acid Sequences)', 'rich_shannon_strip.jpg')
    strip_plot(ov_df, xslice, 'richness_vgene', 
           'Richness (TCR V-gene)', 'rich_vgene_strip.jpg')
    strip_plot(ov_df, xslice, 'richness_jgene', 
           'Richness (TCR J-gene)', 'rich_jgene_strip.jpg')
    strip_plot(ov_df, xslice, 'richness_vdjgene', 
           'Richness (Unique V, D, J-gene Combinations)', 
           'rich_vdjgene_strip.jpg')