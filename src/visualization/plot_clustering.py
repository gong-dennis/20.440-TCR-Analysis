import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
import pickle as pkl
from src.visualization.plotting import strip_plot, ov_strip_plot_stats

def survival_vs_deepcat(ov_df):
    ov_slice = ov_df[['sample_name', 'survival_months', 'group_label']]
    prefix = 'data/processed/shared_seqs_deepcat/'

    db_matches = _get_mults(prefix)
    ov_slice['db_matches'] = db_matches
    
    ax = sns.scatterplot(data=ov_slice, x='db_matches', y='survival_months', 
                         hue='group_label')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.xscale('log')
    plt.xlabel('Log of Number of Cancer-Associated TCRs')
    plt.ylabel('Months of Survival')
    plot_prefix = 'fig/main_fig/cluster_plots/'
    plt.savefig(plot_prefix + 'survival_vs_matches.jpg', 
                bbox_inches='tight')

    ov_strip_plot_stats(ov_slice, 'group_label', 'survival_months', 
               'Months of Survival', plot_prefix + 'survivals_nocluster.jpg')
    ov_strip_plot_stats(ov_slice, 'group_label', 'db_matches', 
               'Number of Cancer-Associated TCRs', 
               plot_prefix + 'matches_nocluster.jpg')
    

def survival_vs_deepcat_clustering(ov_df):
    ov_slice = ov_df[['sample_name', 'survival_months', 'group_label']]
    prefix = 'data/processed/shared_seqs_deepcat/'

    db_matches = _get_mults_with_clustering(prefix)
    ov_slice['db_matches'] = db_matches
    ax = sns.scatterplot(data=ov_slice, x='db_matches', y='survival_months', 
                         hue='group_label')
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.xscale('log')
    plt.xlabel('Log of Number of Cancer-Associated TCRs')
    plt.ylabel('Months of Survival')
    plot_prefix = 'fig/main_fig/cluster_plots/'
    plt.savefig(os.path.join('fig', 'main_fig', 'cluster_plots', 
                             'survival_vs_matches_clustering.jpg'), 
                             bbox_inches='tight')
    
    ov_strip_plot_stats(ov_slice, 'group_label', 'survival_months', 
               'Months of Survival', plot_prefix + 'survivals_cluster.jpg')
    ov_strip_plot_stats(ov_slice, 'group_label', 'db_matches', 
               'Number of Cancer-Associated TCRs', 
               plot_prefix + 'matches_cluster.jpg')


def _get_mults_with_clustering(cancer_prefix):
    mults_prefix = 'data/processed/multiplicity_dicts/'
    cluster_prefix = 'data/processed/clustcr_labels/'
    fnames = sorted(os.listdir(cancer_prefix))
    mults = []
    for f in fnames:
        f_handle = os.path.splitext(f)[0]
        mults_dict = pkl.load(open(mults_prefix + f_handle + '.p', 'rb'))
        cancer_tcrs = pd.read_csv(cancer_prefix + f).iloc[:,0].to_list()
        clusters = pd.read_csv(cluster_prefix + f)

        cum_mult = 0
        for seq in cancer_tcrs:
            if seq in clusters['junction_aa'].values:
                slice = clusters.loc[clusters['junction_aa'] == seq]
                cluster_num = slice['cluster'].iloc[0]
                cluster_members = clusters.loc[
                    clusters['cluster'] == cluster_num]['junction_aa'].to_list()
                for c in cluster_members:
                    cum_mult += mults_dict[c]
            else:
                cum_mult += mults_dict[seq]
        mults.append(cum_mult)
    return mults


def _get_mults(cancer_prefix):
    mults_prefix = 'data/processed/multiplicity_dicts/'
    fnames = sorted(os.listdir(cancer_prefix))
    mults = []
    for f in fnames:
        f_handle = os.path.splitext(f)[0]
        mults_dict = pkl.load(open(mults_prefix + f_handle + '.p', 'rb'))
        cancer_tcrs = pd.read_csv(cancer_prefix + f).iloc[:,0].to_list()
        cum_mult = 0
        for aa in cancer_tcrs:
            cum_mult += mults_dict[aa]
        mults.append(cum_mult)
    return mults