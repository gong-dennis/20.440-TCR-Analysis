import seaborn as sns
import matplotlib.pyplot as plt
import os

def violin(ov_df, xslice, yslice, xlabel, ylabel, filename, show=False):
    plt.figure()
    sns.violinplot(data=ov_df, x=xslice, y=yslice)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(os.path.join('fig', 'main_fig', filename), bbox_inches='tight')
    
    if show: plt.show()


def strip_plot(ov_df, xslice, yslice, ylabel, filename, show=False):
    plt.figure()
    sns.stripplot(data=ov_df, x=xslice, y=yslice)
    plt.xlabel('')
    plt.ylabel(ylabel)

    plt.savefig(os.path.join('fig', 'main_fig', 'overview_strip_plots', 
                             filename), bbox_inches='tight')
    
    if show: plt.show()