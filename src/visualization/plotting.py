import seaborn as sns
import matplotlib.pyplot as plt
import os
from statannotations.Annotator import Annotator
from scipy.stats import mannwhitneyu

def violin(ov_df, xslice, yslice, xlabel, ylabel, filename, show=False):
    plt.figure()
    sns.violinplot(data=ov_df, x=xslice, y=yslice)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(os.path.join('fig', 'main_fig', filename), bbox_inches='tight')
    
    if show: plt.show()


def strip_plot(ov_df, xslice, yslice, ylabel, filepath, show=False):
    plt.figure()
    sns.stripplot(data=ov_df, x=xslice, y=yslice)
    plt.xlabel('')
    plt.ylabel(ylabel)

    plt.savefig(filepath, bbox_inches='tight')
    
    if show: plt.show()


def ov_strip_plot_stats(ov_df, xslice, yslice, ylabel, filepath, show=False, 
                        y_log=False):
    plt.clf()
    plt.figure(figsize=(4,6))
    
    plotting_parameters = {
        'data': ov_df,
        'x': xslice,
        'y': yslice,
    }

    pairs = [('No NACT', 'Short Interval'),
            ('No NACT', 'Long Interval'),
            ('Short Interval', 'Long Interval')]

    with sns.plotting_context('notebook'):
        # Plot with seaborn
        ax = sns.stripplot(**plotting_parameters)
        if y_log: ax.set(yscale="log")

        no_nact = ov_df.loc[ov_df['group_label']=='No NACT', yslice]
        short_int = ov_df.loc[ov_df['group_label']=='Short Interval', yslice]
        long_int = ov_df.loc[ov_df['group_label']=='Long Interval', yslice]
        pvalues = [
            mannwhitneyu(no_nact, short_int, alternative="two-sided").pvalue,
            mannwhitneyu(no_nact, long_int, alternative="two-sided").pvalue,
            mannwhitneyu(short_int, long_int, alternative="two-sided").pvalue
        ]

        # Transform each p-value to "p=" in scientific notation
        formatted_pvalues = _annotate_p_vals(pvalues, alpha=0.05)

        indices = [i for i, x in enumerate(formatted_pvalues) if x == 'ns']

        if not len(indices) == 3:
            for index in sorted(indices, reverse=True):
                del pairs[index]
                del formatted_pvalues[index]

        # Add annotations
        annotator = Annotator(ax, pairs, **plotting_parameters)
        annotator.set_custom_annotations(formatted_pvalues)
        annotator.annotate()

        plt.xlabel('')
        plt.ylabel(ylabel)

        if show: plt.show()
        plt.savefig(filepath, bbox_inches='tight')


def _annotate_p_vals(pvalues, alpha=0.05):
    annotations = []
    for p in pvalues:
        if p < 0.0001:
            annotations.append('****')
        elif p < 0.001:
            annotations.append('***')
        elif p < 0.01:
            annotations.append('**')
        elif p < 0.05:
            annotations.append('*')
        else:
            annotations.append('ns')
    
    return annotations





        