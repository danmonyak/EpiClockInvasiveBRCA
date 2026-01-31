"""
util.py
=======
Author - Daniel Monyak
6-7-24
=======

Master module for this project
Functions and variables used by all scripts and notebooks

"""

import numpy as np
import pandas as pd
import os
import sys
from math import ceil, isnan
from itertools import product, accumulate
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress, ranksums
import json

############################################################
"""
Find the absolute path to the repository directory (repo_dir) whenever this module is imported
repo_dir is needed to read src/consts.json
"""
################################################
exception_msg = '-'*120
exception_msg += '\nutil.py module can only imported inside of the EpiClockInvasiveBRCA directory (or a subdirectory)'
exception_msg += '\nIf running a .py script, navigate to the subdirectory and run "python <script.py>"'
exception_msg += '\nIf running a .ipynb jupyter notebook, ensure that the jupyter server was started inside the EpiClockInvasiveBRCA directory.'
exception_msg += '\n' + '-'*139

subdir_list = os.getcwd().split(os.sep)
while True:
    try:
        if subdir_list[-1] == 'EpiClockInvasiveBRCA':
            break
    except IndexError:
        raise Exception(exception_msg) from None
    
    subdir_list.pop()


repo_dir = os.path.join(os.sep, *subdir_list)
################################################
# Load the variables in consts.json into a dictionary
consts = json.loads(''.join(open(os.path.join(repo_dir, 'src', 'consts.json'), 'r').readlines()))
consts['repo_dir'] = repo_dir
try:
    config = json.loads(''.join(open(os.path.join(repo_dir, 'config.json'), 'r').readlines()))
except FileNotFoundError:
    sys.exit('Please create config.json file in main directory of repository...')
consts.update(config)
################################################


################################################################
##########           Small helper functions           ##########
################################################################

def splitAndJoin(x, n):
    """
    Return the first n fields of a string (separated by dashes)
    """
    return '-'.join(x.split('-')[:n])
def sampleToPatientID(x):
    return splitAndJoin(x, 3)
def getSampleID(x):
    return splitAndJoin(x, 4)

# Element-wise boolean test - is the element None or NaN
isNaVec = np.vectorize(lambda x:(x is None) or ((type(x) is not str) and isnan(x)))

def combineFilters(filters):
    """
    Combine a list of boolean arrays using the bitwise AND operation
    """
    return list(accumulate(filters, lambda x,y:x&y))[-1]

################################################################
##########        Statistical helper functions        ##########
################################################################

def pearsonCorrelation(ser1, ser2, get_n_used=False):
    """
    Compute Pearson correlation between two Pandas Series that represent two vectors
    
    Parameters
    ----------
    ser1, ser2 : Pandas Series objects of numerical dtype
        Two vectors to test correlation of
    get_n_used : boolean
        True iff user wants to also return "n"
    
    Returns
    -------
    res or (res, use_mask)
    res : LinregressResult instance
        Pearson correlation results (output of linregress)
    n : numerical
        Number of pairs of values used for the correlation
    """
    use_mask = ~(ser1.isna() | ser2.isna())
    res = linregress(ser1[use_mask], ser2[use_mask])
    if get_n_used:
        n = use_mask.sum()
        return res, n
    else:
        return res

def getCorrelation(sample_annotations, var_x, var_y, use_samples=None, get_n_used=False):
    """
    Compute Pearson correlation between two columns of a DataFrame
    
    Parameters
    ----------
    sample_annotations : Pandas Dataframe
        has columns var_x and var_y
    var_x : str
        numerical variable x
    var_y : str
        numerical variable y
    use_samples : list or ndarray of strs
        samples to restrict correlation to
    get_n_used : boolean
        True iff user wants to also return "n"
    
    Returns
    -------
    res or (res, use_mask)
    res : LinregressResult instance
        Pearson correlation results (output of linregress)
    n : numerical
        Number of pairs of values used for the correlation
    """
    mask = ~sample_annotations[var_x].isna()
#     if only_pure:
#         mask &= sample_annotations['pure']
    if use_samples is not None:
        mask &= sample_annotations.index.isin(use_samples)
    df = sample_annotations.loc[mask]
    ser1 = df[var_y]
    ser2 = df[var_x]
    return pearsonCorrelation(ser1, ser2, get_n_used)

def wilcoxonRankSums(ser1, ser2):
    """
    Compute the Wilcoxon rank-sum statistic for two samples held in Pandas Series
    
    Parameters
    ----------
    ser1, ser2 : Pandas Series objects of numerical dtype
    
    Returns
    -------
    Wilcoxon rank-sum results (output of ranksums)

    """
    return ranksums(ser1.dropna(), ser2.dropna())

def getWilcoxonPvalueTable(sample_annotations, var_cat, var_y, use_groups=None):
    """
    Compute the Wlcoxon rank-sum statistic pvalue between all combinations of possible
        groups of a categorical variable in a DataFrame
    
    Parameters
    ----------
    sample_annotations : Pandas Dataframe
        has columns var_cat and var_y
    var_cat : str
        categorical variable to stratify by
    var_y : str
        numerical variable to measure
    use_groups : list of strs
        specific values of var_cat to compare
        if None, use all unique values of var_cat
    
    Returns
    -------
    pvalue_df : Pandas DataFrame
        Wilcoxon rank-sum p-values for each pair of values of var_cat

    """
    if use_groups is None:
        use_groups = sample_annotations[var_cat].unique()
        
    pvalue_df = pd.DataFrame(index=use_groups, columns=use_groups, data=float('nan'))

    for i in range(len(use_groups)):
        group_i = use_groups[i]
        for j in range(i+1, len(use_groups)):
            group_j = use_groups[j]
            ser_i = sample_annotations.loc[sample_annotations[var_cat] == group_i, var_y]
            ser_j = sample_annotations.loc[sample_annotations[var_cat] == group_j, var_y]
            res = wilcoxonRankSums(ser_i, ser_j)
            pvalue_df.loc[group_i, group_j] = res.pvalue

    return pvalue_df

#############################################################
#########         LUMP purity calculation         ###########
#############################################################

def getLUMP_values(beta_values):
    """
    Compute the LUMP purty value of each sample from their LUMP site beta values
    min(1, average beta value of LUMP sites divided by 0.85)
    
    Parameters
    ----------
    beta_values : Pandas Dataframe
        rows are CpGs, columns are samples
    
    Returns
    -------
    pvalue_df : Pandas DataFrame
        Wilcoxon rank-sum p-values for each pair of values of var_cat

    """
    
    # Import 44 LUMP CpG site names
    lump_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], 'data', 'lump-CpGs-44.txt'), dtype=str)
    
    included_lump = np.intersect1d(lump_CpGs, beta_values.index)
    if included_lump.shape[0] != lump_CpGs.shape[0]:
        print(f'Only {included_lump.shape[0]} LUMP sites were available...')
    raw_LUMPs = (beta_values.loc[included_lump].mean(axis=0) / 0.85).to_frame()
    raw_LUMPs[1] = 1
    return raw_LUMPs.min(axis=1)


###############################################################
##########         Plotting helper functions         ##########
###############################################################

def saveBoxPlotNew(sample_annotations, var_cat, var_y='c_beta', restrict=True, use_groups=None,
                outdir='.', outfile=True, title=False, custom_title=None, xlabel=None, ylabel=None, palette=None,
                plot_ymax_mult=0.25, signif_bar_heights=0.03,
                   signif_fontsize=14, ylim=None,
                   figsize=(10, 10), labelfontsize=20, ticksfontsize=10, linewidth=1, fliersize=1, sf=1):
    """
    Create box plot from a sample annotations DataFrame
    Plot some numerical variable and stratfy by some categorical variable
    Can include significance bars (Wilcoxon rank-sum)
    
    Parameters
    ----------
    sample_annotations : Pandas Dataframe
        has columns var_cat and var_y
    var_cat : str
        categorical variable to stratify by
    var_y : str
        numerical variable to measure
    restrict : boolean
        True iff we should restrict to samples with in_analysis_dataset==True
    use_groups : list of strs
        specific values of var_cat to compare
        if None, use all unique values of var_cat
    outdir : str
        output directory of figure file generated
    outfile : boolean
        True iff we want to save the figure to a file
    title : boolean
        use title
    custom_title : str
        title to use other than default (sample_annotations.name)
    xlabel : str
        custom x-label
    ylabel : str
        custom y-label
    palette : palette name, list, or dict
        color palette
    plot_ymax_mult : float
        proportionally changes y top limit
    signif_bar_heights : float
        height of vertical part of significance bars
        set to None if you don't want significance bars
    signif_fontsize : float
        fontsize of signficances * symbols
    ylim : tuple of floats
        y-axis limits that is used if not None and signif_bar_heights is None
    figsize : tuple of floats
        figure size
    labelfontsize : float
        fontsize of x-label, y-label, x-ticks (categories), and title
    ticksfontsize : float
        fontsize of y-ticks
    linewidth : float
        linewidth argument passed to sns.boxplot
    fliersize : float
        fliersize argument passed to sns.boxplot
    sf : float
        scale factor
        change to alter the size of a figure - scales everything proportionally
    
    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes instance used for plot

    """
    
    ########################
    ##### Select data
    ########################
    
    # Select samples wth var_cat in use_groups
    if use_groups is None:
        use_groups = sample_annotations[var_cat].unique()
        use_groups = np.sort(use_groups[~isNaVec(use_groups)])
    
    # Boolean mask of samples to use
    use_samples_mask = sample_annotations[var_cat].isin(use_groups)

    # Add to mask and set outfile name
    if restrict:
        use_samples_mask &= sample_annotations['in_analysis_dataset']
        if outfile:
            outfile_name = f'{sample_annotations.name}-{var_cat}-{var_y}-restrict.pdf'
    elif outfile:
        outfile_name = f'{sample_annotations.name}-{var_cat}-{var_y}.pdf'

    # Exclude samples with missing data
    use_samples_mask &= ~sample_annotations[var_y].isna()
    
    # Define final DataFrame for plotting
    plot_data = sample_annotations.loc[use_samples_mask, [var_cat, var_y]]
    
    # Sanity check
    assert plot_data.shape[0] == plot_data.dropna().shape[0]
    
    ########################
    ##### Create plot
    ########################
    
    fig, ax = plt.subplots(figsize=np.array(figsize) * sf)
    sns.boxplot(ax=ax, data=plot_data, x=var_cat, y=var_y,
                order=use_groups, palette=palette,
                linewidth=linewidth * sf, fliersize=fliersize * sf)

    ######################################################
    ##### Customize plot labels, axes, ticks, ticklabels
    ######################################################
    
    if var_y == 'c_beta':
        ax.set_ylabel('$c_β$', fontsize=labelfontsize * sf)
    elif var_y == 'c_beta_adj1':
        ax.set_ylabel('$c_β^α$', fontsize=labelfontsize * sf)
    elif ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=labelfontsize * sf)

    ax.tick_params(axis='x', labelsize=labelfontsize * sf, width=sf, length=8 * sf)
    ax.tick_params(axis='y', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)

    # Add (n = ...) under each x-tick (category)
    ax.set_xticks(ax.get_xticks(),
                  [str(group) + f'\n(n = {(~(plot_data[var_y].isna()) & (plot_data[var_cat] == group)).sum()})' for group in use_groups])
    
    # Set title
    # Default = sample_annotations.name
    # By default, xlabel is erased from plot
    if title:
#         if label is None:
        ax.set_title(sample_annotations.name, fontsize=labelfontsize * sf)
#         else:
#             ax.set_title(f'{label} ({sample_annotations.name})', fontsize=labelfontsize * sf)
        ax.set_xlabel('', fontsize=labelfontsize * sf)
    elif custom_title is not None:
        ax.set_title(custom_title, fontsize=labelfontsize * sf)
    elif xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=labelfontsize * sf)

    ######################################################
    ##### Add significance bars
    ######################################################
    
    if signif_bar_heights is not None:
        max_y = plot_data[var_y].max()
        min_y = plot_data[var_y].min()

        # Create list of indexes
        idxs_1 = list(range(1, len(use_groups) + 1))
        
        # Create list of combinations of indexes
        pair_list = [(idxs_1[x], idxs_1[x + y]) for y in reversed(idxs_1) for x in range((len(idxs_1) - y))]

        # Hold combinations with a signifcant p-value
        significant_pairs = []
        for combo in pair_list:

            ser1 = plot_data.loc[plot_data[var_cat] == use_groups[combo[0]-1], var_y]
            ser2 = plot_data.loc[plot_data[var_cat] == use_groups[combo[1]-1], var_y]

            pvalue = wilcoxonRankSums(ser1, ser2).pvalue
            if pvalue < 0.05:
                significant_pairs.append([combo, pvalue])

        plot_ymax = plot_ymax_mult * 0.8 * max(1, len(significant_pairs)) * (max_y - min_y) + max_y
        ax.set_ylim(top=plot_ymax)

        # Get the y-axis limits
        bottom, top = ax.get_ylim()
        # top = plot_ymax
        y_range = top - bottom

        # Significance bars
        for i, combo in enumerate(significant_pairs):
            # Columns corresponding to the datasets of interest
            x1 = combo[0][0] - 1
            x2 = combo[0][1] - 1
            
            level = len(significant_pairs) - i
            
            # Plot the bar
            bar_height = max_y + (signif_bar_heights * level)

            bar_tips = bar_height - (y_range * 0.02)
            ax.plot(
                [x1, x1, x2, x2],
                [bar_tips, bar_height, bar_height, bar_tips], lw=linewidth * sf, c='k'
            )
            # Significance level
            pvalue = combo[1]
            if pvalue < 0.001:
                sig_symbol = '***'
            elif pvalue < 0.01:
                sig_symbol = '**'
            elif pvalue < 0.05:
                sig_symbol = '*'
            text_height = bar_height
            ax.text((x1 + x2) * 0.5, text_height, sig_symbol,
                    ha='center', va='bottom', c='k', fontsize=signif_fontsize * sf)
    elif ylim is not None:
        ax.set_ylim(ylim)
    
    # Save plot
    if outfile:
        fig.savefig(os.path.join(outdir, outfile_name), format='pdf', pad_inches=0.1)

    return ax

def saveCorrelationPlot(sample_annotations, var_y, var_x='c_beta', restrict=True, use_samples=None,
                        outdir='.', outfile=True, text_x=0, text_y=0,
                        scatter_kws={}, line_kws={}, ylabel=None, bbox_dict=None, color='blue',
                       figsize=(10, 10), labelfontsize=20, ticksfontsize=10, s=1, sf=1):
    """
    Create scatter plot with a line of best fit from a sample annotations DataFrame
    Plot two numerical variables against each other
    
    Parameters
    ----------
    sample_annotations : Pandas Dataframe
    var_y : str
        numerical variable to plot on y-axis
    var_x : str
        numerical variable to plot on x-axis
    restrict : boolean
        True iff we should restrict to samples with in_analysis_dataset==True
    use_samples : list or ndarray of strs
        samples to restrict correlation to
    outdir : str
        output directory of figure file generated
    outfile : boolean
        True iff we want to save the figure to a file
    text_x : float
        x coordinate of "R = ..." text
    text_y : float
        y coordinate of "R = ..." text
    scatter_kws : dict
        scatter plot keyword arguments -- see sns.regplot documentation
    line_kws : dict
        line plot keyword arguments -- see sns.regplot documentation
    ylabel : str
        custom y-label
        default is var_y
    color : matplotlib color
        color of all plot elements
    figsize : tuple of floats
        figure size
    labelfontsize : float
        fontsize of x-label, y-label, and title
    ticksfontsize : float
        fontsize of x and y-ticks
    s : float
        size of scatter points
    sf : float
        scale factor
        change to alter the size of a figure - scales everything proportionally
        
    Notes
    -----
    sample_annotations.name must be set
    i.e. sample_annotations.name = 'dataset'

    """
    
    ########################
    ##### Select data
    ########################
    
    # Boolean mask of samples to use
    use_samples_mask = np.ones(shape=sample_annotations.shape[0], dtype=bool)
    
    # Add to mask and set outfile name
    if restrict:
        use_samples_mask &= sample_annotations['in_analysis_dataset']
        outfile_name = f'{sample_annotations.name}-{var_x}-{var_y}-pure.pdf'
    else:
        outfile_name = f'{sample_annotations.name}-{var_x}-{var_y}.pdf'
    
    # Set use_samples if it is None
    if use_samples is None:
        use_samples = sample_annotations.index[use_samples_mask]
    else:
        use_samples = np.intersect1d(use_samples, sample_annotations.index[use_samples_mask])
    
    # Define final DataFrame for plotting
    plot_data = sample_annotations.loc[use_samples, [var_x, var_y]].dropna()
    
    # Calculate Pearson correlation between variables
    res = getCorrelation(plot_data, var_x=var_x, var_y=var_y)
    
    ########################
    ##### Create plot
    ########################
    
    fig, ax = plt.subplots(figsize=np.array(figsize) * sf)
    scatter_kws['s'] = s * sf**2
    sns.regplot(ax=ax, data=plot_data, x=var_x, y=var_y, scatter_kws=scatter_kws, color=color,
               line_kws=line_kws)

    ######################################################
    ##### Customize plot labels, axes, ticks, ticklabels
    ######################################################

    # Y-label by default is var_y, split by _ and capitalized
    if ylabel is None:
        ax.set_ylabel(" ".join(var_y.split("_")).capitalize(), fontsize=labelfontsize * sf)
    else:
        ax.set_ylabel(ylabel, fontsize=labelfontsize * sf)

    ax.set_title(sample_annotations.name, fontsize=labelfontsize * sf)
    ax.tick_params(axis='both', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)
    
    # Put (R = ...) on the plot
    ax.text(text_x, text_y, f'R = {res.rvalue:.2f}',
                        ha="center", va="bottom",
                        fontfamily='sans-serif', fontsize=0.8 * labelfontsize * sf, bbox=bbox_dict)
    
    if var_x == 'c_beta':
        ax.set_xlabel('$c_β$', fontsize=labelfontsize * sf)
    elif var_x == 'c_beta_adj1':
        ax.set_xlabel('$c_β^α$', fontsize=labelfontsize * sf)
        
    # Save plot
    if outfile:
        fig.savefig(os.path.join(outdir, outfile_name), format='pdf', pad_inches=0.1)


def plotTumorWise(beta_values, CpG_list=None, sample_list=None, n_samps=30, ncols=3, suptitle='random pick of samples',
                  title_formats=None, xlabel='Beta', random_seed=None, outfile=False, outfile_name=None,
                  outdir='images', choose_random=True, color='blue', ylim=None, bins='auto', figsize=None,
                  text_fontsize=None, ticksfontsize=None, opacity=None, sf=1, tight_layout_pad=1, kde=False,
                  suptitle_y=0.99, suptitle_fontsize=20):
    """
    Create a panel of tumor-level histograms of beta values
    Can specify a list of CpGs to plot and/or a list of samples to plot
    Can randomly pick samples to plot
    
    Parameters
    ----------
    beta_values : Pandas Dataframe
        rows are CpGs, columns are samples
    CpG_list : list or ndarray of strs
        CpG sites to plot for each sample
    sample_list : list or ndarray of strs
        samples to plot
    n_samps : int
        number of samples to plot
    ncols : int
        number of columns
    suptitle : str
        super title of plot
    title_formats : str
        string with one "{}" instance to place sample names in
        custom format of title of each hstogram
    xlabel : str
        x-axis title
    random_seed : numerical
        random_seed to use when randomly picking samples
    outfile : boolean
        True iff we want to save the figure to a file
    outfile_name : str
        name of outfile
    outdir : str
        path to output directory
    choose_random : boolean
        True iff we want to randomly pick samples
    color : matplotlib color or list of matplotlib colors
        color of histogram
        if list, use colors in order for the first len(color) plots,
        then use last color in list if it needs more
    ylim : tuple of floats
        y-axis limits that is used if not None
    bins : str, number, vector, or a pair of such values
        bins argument to sns.histplot
    figsize : tuple of floats
        figure size
    text_fontsize=float
        fontsize of x-label, y-label, and title
    ticksfontsize : float
        fontsize of x and y-ticks
    opacity : float
        opacity level to pass to alpha argument in sns.histplot
    sf : float
        scale factor
        change to alter the size of a figure - scales everything proportionally
    tight_layout_pad : float
        tight_layout_pad argument to pass to fig.tight_layout
    kde : boolean
        True iff it should plot a kde curve
    suptitle_y : float
        y parameter in fig.suptitle
        relative height of suptitle
    suptitle_fontsize : float
        fontsize of suptitle
    """
    
    ########################
    ##### Select data
    ########################
    
    # Select samples from sample_list or columns of beta_values
    if sample_list is None:
        sample_list = beta_values.columns.values
    
    n_samps = min(n_samps, len(sample_list))
    nrows = ceil(n_samps / ncols)
    
    # Use a list of CpG sites or all available in beta_values
    if CpG_list is None:
        CpG_list = beta_values.index.values
    
    # Choose a random set of samples
    if choose_random:
        np.random.seed(random_seed)
        samples_randSamp = np.random.choice(sample_list, n_samps, replace=False)
    else:
        samples_randSamp = sample_list[:n_samps]
    
    ########################
    ##### Create plot
    ########################
    
    if figsize is None:
        fig, axes_arr = plt.subplots(nrows, ncols, figsize=(20, 3 + 3*nrows))
    else:
        fig, axes_arr = plt.subplots(nrows, ncols, figsize=np.array(figsize) * sf)
    fig.suptitle(suptitle, y=suptitle_y, fontsize=suptitle_fontsize, fontweight='bold')
    fig.tight_layout(pad=tight_layout_pad)
    
    # Iterate and plot all samples
    for i, samp in enumerate(samples_randSamp):
        col = i % ncols
        if nrows > 1:
            row = i // ncols
            ax = axes_arr[row, col]
        else:
            ax = axes_arr[col]
        
        if type(color) is list:
            cur_color = color[min(i, n_samps-1)]
        else:
            cur_color = color

        plot = sns.histplot(ax=ax, data=beta_values.loc[CpG_list, samp],
                            stat='proportion', binrange=(0, 1),
                           color=cur_color, bins=bins, alpha=opacity, kde=kde)
        
        # Set title of local histogram
        if title_formats is None:
            title = samp
        else:
            title = title_formats[i].format(samp)
        
        # Customize plot labels, axes, ticks, ticklabels
        ax.set_title(title, fontsize=text_fontsize * sf)
        ax.set_xlabel(xlabel, fontsize=text_fontsize * sf)
        if col == 0:
            ax.set_ylabel('Proportion', fontsize=text_fontsize * sf)
        else:
            ax.set_ylabel('')
        ax.tick_params(axis='both', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)
        
        if ylim is not None:
            ax.set_ylim(ylim)
    
    if not outfile:
        fig.show()
    elif outfile_name is None:
        print('Provide a file name...')
    else:      # Save plot
        fig.savefig(os.path.join(outdir, outfile_name), format='pdf', pad_inches=0.1)    



def saveFigureData(figure_data, figure_name):
    """
    Save a dataframe in the Figure_data_dir directory, WITHOUT the index, to the "Figure_data_dir" directory
    
    Parameters
    ----------
    figure_data : Pandas Dataframe
        data to save
    figure_name : str
        prefix of filename
    """
    figure_data.to_csv(os.path.join(consts['Figure_data_dir'], f'{figure_name}.tsv'), sep='\t',
                      index=False)
