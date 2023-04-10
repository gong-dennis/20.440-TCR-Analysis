from src.visualization.plotting import strip_plot, ov_strip_plot_stats

def create_stripplots(ov_df):
    ## Plotting
    prefix = 'fig/main_fig/overview_strip_plots/'
    xslice = 'group_label'
    strip_plot(ov_df, xslice, 'fraction_productive', 
           'Fraction Productive', prefix + 'fraction_productive_strip.jpg')
    strip_plot(ov_df, xslice, 'rel_prod_rearr', 
           'Relative Productive Rearrangements', prefix + 'rel_prod_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'rel_total_rearr', 
           'Relative Total Rearrangements', prefix + 'rel_total_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'frac_prod_rearr', 
           'Fraction of Productive Rearrangements', prefix + 'frac_prod_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'productive_simpson_clonality', 
           'Productive Simpson Clonality', prefix + 'prod_simps_clon_rearr_strip.jpg')
    strip_plot(ov_df, xslice, 'max_productive_frequency', 
           'Max Productive Frequency', prefix + 'max_prod_freq_strip.jpg')
    
    strip_plot(ov_df, xslice, 'aa_shannon', 
           'Shannon Index (from Amino Acid Sequences)', prefix + 'aa_shannon_strip.jpg')
    strip_plot(ov_df, xslice, 'mean_cdr3', 
           'Mean CDR3 Length', prefix + 'mean_cdr3_strip.jpg')
    strip_plot(ov_df, xslice, 'shannon_vgene', 
           'Shannon Index (TCR V-gene)', prefix + 'shannon_vgene_strip.jpg')
    strip_plot(ov_df, xslice, 'shannon_jgene', 
           'Shannon Index (TCR J-gene)', prefix + 'shannon_jgene_strip.jpg')
    strip_plot(ov_df, xslice, 'shannon_vdjgene', 
           'Shannon Index (Unique V, D, J-gene Combinations)', 
           prefix + 'shannon_vdjgene_strip.jpg')
    
    strip_plot(ov_df, xslice, 'aa_gini', 
           'Gini Impurity (from Amino Acid Sequences)', prefix + 'aa_gini_strip.jpg')
    strip_plot(ov_df, xslice, 'gini_vgene', 
           'Gini Impurity (TCR V-gene)', prefix + 'gini_vgene_strip.jpg')
    strip_plot(ov_df, xslice, 'gini_jgene', 
           'Gini Impurity (TCR J-gene)', prefix + 'gini_jgene_strip.jpg')
    strip_plot(ov_df, xslice, 'gini_vdjgene', 
           'Gini Impurity (Unique V, D, J-gene Combinations)', 
           prefix + 'gini_vdjgene_strip.jpg')
    
    strip_plot(ov_df, xslice, 'aa_richness', 
           'Richness (from Amino Acid Sequences)', prefix + 'rich_shannon_strip.jpg')
    strip_plot(ov_df, xslice, 'richness_vgene', 
           'Richness (TCR V-gene)', prefix + 'rich_vgene_strip.jpg')
    strip_plot(ov_df, xslice, 'richness_jgene', 
           'Richness (TCR J-gene)', prefix + 'rich_jgene_strip.jpg')
    strip_plot(ov_df, xslice, 'richness_vdjgene', 
           'Richness (Unique V, D, J-gene Combinations)', 
           prefix + 'rich_vdjgene_strip.jpg')
    

    # ------
    ov_strip_plot_stats(ov_df, xslice, 'aa_shannon', 
                        'Shannon Index (from Amino Acid Sequences)', 
                        prefix + '_stats_aa_shannon_strip.jpg')
    ov_strip_plot_stats(ov_df, xslice, 'shannon_vdjgene', 
                        'Shannon Index (Unique V, D, J-gene Combinations)', 
                        prefix + '_stats_shannon_vdjgene_strip.jpg')
    ov_strip_plot_stats(ov_df, xslice, 'aa_richness', 
                        'Richness (from Amino Acid Sequences)', 
                        prefix + '_stats_rich_shannon_strip.jpg')
    ov_strip_plot_stats(ov_df, xslice, 'richness_vdjgene', 
                        'Richness (Unique V, D, J-gene Combinations)', 
                        prefix + '_stats_rich_vdjgene_strip.jpg')

    