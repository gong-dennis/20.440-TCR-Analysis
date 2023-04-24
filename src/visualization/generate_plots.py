from os.path import join
import pandas as pd

from src.visualization.plotting import strip_plot, ov_strip_plot_stats

def plot_pop_metrics():
    ov_metric_path = 'data/processed/overview_with_metrics.csv'
    df = pd.read_csv(ov_metric_path, sep="\t")

    prefix = 'fig/main_fig/pop_metrics_strip'
    xslice = 'group_label'
    strip_plot(df, xslice, 'richness_aa', 'Richness (AA Sequences)', 
           join(prefix, 'richness_aa.jpg'))
    strip_plot(df, xslice, 'richness_vdjgene', 'Richness (Unique VDJ)', 
           join(prefix, 'richness_vdjgene.jpg'))
    strip_plot(df, xslice, 'gini_aa', 'Gini Impurity (AA Sequences)', 
           join(prefix, 'gini_aa.jpg'))
    strip_plot(df, xslice, 'gini_vdjgene', 'Gini Impurity (Unique VDJ)', 
           join(prefix, 'gini_vdjgene.jpg'))
    strip_plot(df, xslice, 'shannon_aa', 'Shannon Index (AA Sequences)', 
           join(prefix, 'shannon_aa.jpg'))
    strip_plot(df, xslice, 'shannon_vdjgene', 'Shannon Index (Unique VDJ)', 
           join(prefix, 'shannon_vdjgene.jpg'))
    
    ov_strip_plot_stats(df, xslice, 'richness_aa', 'Richness (AA Sequences)', 
           join(prefix, 'stats_richness_aa.jpg'))
    ov_strip_plot_stats(df, xslice, 'richness_vdjgene', 'Richness (Unique VDJ)', 
           join(prefix, 'stats_richness_vdjgene.jpg'))
    ov_strip_plot_stats(df, xslice, 'gini_aa', 'Gini Impurity (AA Sequences)', 
           join(prefix, 'stats_gini_aa.jpg'))
    ov_strip_plot_stats(df, xslice, 'gini_vdjgene', 'Gini Impurity (Unique VDJ)', 
           join(prefix, 'stats_gini_vdjgene.jpg'))
    ov_strip_plot_stats(df, xslice, 'shannon_aa', 'Shannon Index (AA Sequences)', 
           join(prefix, 'stats_shannon_aa.jpg'))
    ov_strip_plot_stats(df, xslice, 'shannon_vdjgene', 'Shannon Index (Unique VDJ)', 
           join(prefix, 'stats_shannon_vdjgene.jpg'))
    

def plot_cancer_assoc_tcrs():
    ov_metric_path = 'data/processed/overview_with_cancer_tcrs.csv'
    df = pd.read_csv(ov_metric_path, sep="\t")

    prefix = 'fig/main_fig/cancer_assoc_tcrs'
    xslice = 'group_label'
    strip_plot(df, xslice, 'num_cancer_tcrs', 'Number of Cancer-Associated TCRs', 
           join(prefix, 'cancer_assoc_tcrs.jpg'))
    ov_strip_plot_stats(df, xslice, 'num_cancer_tcrs', 'Number of Cancer-Associated TCRs', 
                        join(prefix, 'stats_cancer_assoc_tcrs.jpg'), y_log=True)