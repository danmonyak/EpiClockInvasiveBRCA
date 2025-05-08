"""
util_loyferNormals.py
=======
Author - Daniel Monyak
5-20-24
=======

Functions used by Select_fCpGs notebook

but using basal and luminal epithelial normal samples from Loyfer et al.
    
"""

import pandas as pd
import numpy as np
import os
import sys
import EpiClockInvasiveBRCA.src.util as epi_util
# from EpiClockInvasiveBRCA.src.consts import consts
import EpiClockInvasiveBRCA.src.util as epi_util
consts = epi_util.consts

############################################################
"""
Defines acceptable ranges (inclusive) for each value
"""
############################################################
CLOCK_CRITERIA = {
    'basal':{
        'beta_nans':(0, 20),
        'beta_means':(0.4, 0.6)
    },
    'luminal':{
        'beta_nans':(0, 20),
        'beta_means':(0.4, 0.6)
    },
    'tumor':{
        'beta_nans':(0, 20),
        'beta_means':(0.4, 0.6)
    },
    'normals':{
        'beta_nans':(0, 20),
        'beta_means':(0.4, 0.6)
    }
}
############################################################
############################################################

def getNeutralDNACpGs():
    """
    Returns
    -------
    neutral_DNA_CpG_list : list
        List of "neutral" CpGs
        i.e. not associated with any gene or regulatory feature
    
    Notes
    -----
    Checks the manifests of both the 450K and 850K array
    """
    
    chip_annot_dir = os.path.join(consts['TCGA_datadir'], 'chip_annots')
    manifest_450K = pd.read_table(os.path.join(chip_annot_dir, 'chip450_annot_cleaned.txt'), index_col=0, low_memory=False)
    manifest_850K = pd.read_table(os.path.join(chip_annot_dir, 'chip850_annot_cleaned.txt'), index_col=0, low_memory=False)
    
    # Requirements for neutral CpGs
    neutral_filters = [
        manifest_450K['Regulatory_Feature_Group'].isna(),
        manifest_450K['UCSC_RefGene_Name'].isna(),
        manifest_850K['Regulatory_Feature_Group'].isna(),
        manifest_850K['UCSC_RefGene_Name'].isna(),
    ]
    
    # Combine the filters to select neutral CpGs
    neutral_DNA_CpG_list = manifest_450K.index[epi_util.combineFilters(neutral_filters)].values
    
    # Sanity checks
    neutral_manifest_450K = manifest_450K.loc[neutral_DNA_CpG_list]
    neutral_manifest_850K = manifest_850K.loc[neutral_DNA_CpG_list]
    assert neutral_manifest_450K['Regulatory_Feature_Group'].isna().all()
    assert neutral_manifest_450K['UCSC_RefGene_Name'].isna().all()
    assert neutral_manifest_850K['Regulatory_Feature_Group'].isna().all()
    assert neutral_manifest_850K['UCSC_RefGene_Name'].isna().all()

    return neutral_DNA_CpG_list

def getDataDict():
    """
    Return a dictionary that holds data for the the TCGA and normal cohorts
    
    For each cohort, dictionary will hold
        beta_values - # samples x # CpGs
        pureSamples - list of samples after filtering for purity
            this step is already done in TCGA
    
    Returns
    -------
    data : dictionary
        Dictionary with 2 levels
        First level is between 'tumor' and 'normal'
            
    Notes
    -----
    Importing data can take up to 10 minutes once file is downloaded (if on an external file source)
        But should be shorter than this
    """
    
    # Nested dictionary
    data = {'tumor':{}, 'basal':{}, 'luminal':{}, 'normals':{}}
    
    data['tumor']['beta_values'] = pd.read_table(os.path.join(consts['TCGA_datadir'], 'TCGA.methyl.tsv'), index_col=0)
#     data['normal']['beta_values'] = pd.read_table(os.path.join(consts['official_indir'], 'Johnson', 'johnson-beta_values.txt'), index_col=0)
    
    ##### import loyfer normals data
    loyfer_indir = os.path.join(consts['official_indir'], 'Loyfer')
    data['basal']['beta_values'] = pd.read_table(os.path.join(loyfer_indir, 'beta_values-Breast-Basal.tsv'), index_col=0)
    data['luminal']['beta_values'] = pd.read_table(os.path.join(loyfer_indir, 'beta_values-Breast-Luminal.tsv'), index_col=0)
    #####
    
    ### import Johnson et al. breast normals samples
    data['normals']['beta_values'] = pd.read_table(os.path.join(consts['official_indir'], 'Johnson', 'johnson-beta_values.txt'), index_col=0)

    # Check the CpG Sites
    idx_tumor = np.sort(data['tumor']['beta_values'].index.values)
    idx_basal = np.sort(data['basal']['beta_values'].index.values)
    idx_luminal = np.sort(data['luminal']['beta_values'].index.values)
    assert np.all(idx_basal == idx_luminal)

    # Use the same order of CpGs in the index
    data['allCpGs'] = np.intersect1d(idx_tumor, idx_basal)
    data['tumor']['beta_values'] = data['tumor']['beta_values'].loc[data['allCpGs']]
    data['basal']['beta_values'] = data['basal']['beta_values'].loc[data['allCpGs']]
    data['luminal']['beta_values'] = data['luminal']['beta_values'].loc[data['allCpGs']]
    data['normals']['beta_values'] = data['normals']['beta_values'].loc[data['allCpGs']]

    # Clip unnecessary accessions suffixes from sample IDs
    data['tumor']['beta_values'] = data['tumor']['beta_values'].rename(lambda x:'-'.join(x.split('-')[:4]), axis='columns')
    
    # Get pure samples for normals
    # Filter by LUMP value
    data['basal']['purity'] = epi_util.getLUMP_values(data['basal']['beta_values'])
    data['luminal']['purity'] = epi_util.getLUMP_values(data['luminal']['beta_values'])
    data['normals']['purity'] = epi_util.getLUMP_values(data['normals']['beta_values'])
    
    # Don't need purity filtering step
    assert (data['basal']['purity'] > 0.9).all()
    assert (data['luminal']['purity'] > 0.9).all()
    
    data['basal']['pureSamples'] = data['basal']['purity'].index[
        data['basal']['purity'] >= consts['lump_threshold_dict']['Johnson']
    ].values
    data['luminal']['pureSamples'] = data['luminal']['purity'].index[
        data['luminal']['purity'] >= consts['lump_threshold_dict']['Johnson']
    ].values
    
    
    data['normals']['pureSamples'] = data['normals']['purity'].index[
        data['normals']['purity'] >= consts['lump_threshold_dict']['Johnson']
    ].values
    
    # All TCGA samples have already been filtered for purity
    data['tumor']['pureSamples'] = data['tumor']['beta_values'].columns.values
    
    return data

def addMeanStdsNans(data):
    """
    Adds data objects to input dictionary
    
    For the following statistics for each CpG across all samples in the data
        mean
        standard deviation
        # NaN values
    Adds beta_values_SELECTION to the dictionary only if it's not there already
    
    Parameters
    ----------
    data : dict
        Dictionary returned by getDataDict
    """
    
    for cohort in ['tumor', 'basal', 'luminal', 'normals']:
        if 'beta_values_SELECTION' not in data[cohort]:
            data[cohort]['beta_values_SELECTION'] = data[cohort]['beta_values'][data[cohort]['pureSamples']]
        
        data[cohort]['beta_means'] = data[cohort]['beta_values_SELECTION'].mean(axis=1)
        data[cohort]['beta_stds'] = data[cohort]['beta_values_SELECTION'].std(axis=1)
        data[cohort]['beta_nans'] = data[cohort]['beta_values_SELECTION'].isna().sum(axis=1)


def getCpG_list(data, criteria, starting_CpG_list=None, n_select=None, sample_rank='tumor', stat_rank='beta_stds', good_end='higher'):
    """
    Return a set of CpGs that fit certain criteria
    
    Parameters
    ----------
    data : dict
        Dictionary returned by getDataDict
    criteria : dict
        Dictionary that defines acceptable ranges (inclusive) for each data type in each cohort
    starting_CpG_list : list
        CpGs to consider initially
    n_select : numerical
        Number of CpGs to select; if None, return all CpGs that fit criteria
    sample_rank : str
        If n_select is not None, use stat_rank from this cohort to rank the tumors
    stat_rank : str
        If n_select is not None, use this stat from the sample_rank cohort to rank the tumors
    good_end : str
        'higher' or 'lower' - which end of stat_rank should be selected if n_select is not None
    
    Returns
    -------
    CpG_list : ndarray
        Array of names of CpGs that fit criteria
        only include those in starting_CpG_list
        limiting to n_select tumors, ranking using stat_rank in the sample cohort in sample_rank
    """
    
    if 'allCpGs' not in data:
        data['allCpGs'] = data['tumor']['beta_values'].index.values
    if starting_CpG_list is None:
        starting_CpG_list = data['allCpGs']
        
    print(f'Starting with {len(starting_CpG_list)} CpGs')
    
    siteFilters = []
    for sample in criteria.keys():
        samp_criteria = criteria[sample]
        for stat in samp_criteria.keys():
            siteFilters.append(data[sample][stat] >= samp_criteria[stat][0])
            siteFilters.append(data[sample][stat] <= samp_criteria[stat][1])
            
    combined_filter = epi_util.combineFilters(siteFilters)
    criteria_CpGs = np.intersect1d(data['allCpGs'][combined_filter], starting_CpG_list)
    if n_select is None:
        return criteria_CpGs
    
    print(f'Picking {n_select} from with {len(criteria_CpGs)} resulting CpGs')
    
    return getHighestRankedCpGs(data, n_select, sample_rank, stat_rank, good_end, criteria_CpGs)

def getHighestRankedCpGs(data, n_select, sample_rank, stat_rank, good_end, criteria_CpGs):
    """
    Returns the highest ranked CpGs according to some criterion
    
    Parameters
    ----------
    data : dict
        Dictionary returned by getDataDict
    n_select : numerical
        Number of CpGs to select
    sample_rank : str
        Cohort name (e.g. "tumor")
        Use stat_rank from this cohort to rank the tumors
    stat_rank : str
        Statistic name (e.g. "beta_stds")
        Use this stat from the sample_rank cohort to rank the tumors
    good_end : str
        "higher" or "lower" - which end of stat_rank should be selected
    criteria_CpGs : ndarray
        Array of names of highest ranked CpGs, ranking using stat_rank in the sample cohort in sample_rank
        Length n_select
    
    Returns
    -------
    CpG_list : ndarray
        Array of names of CpGs that fit criteria
        only include those in starting_CpG_list
        limiting to n_select tumors, ranking using stat_rank in the sample cohort in sample_rank
    """
    return data[sample_rank][stat_rank].loc[criteria_CpGs].sort_values(ascending=(good_end == 'lower')).index[:n_select].values

def gen_CpG_set(data, neutral_DNA_CpG_list, n_select=500):
    """
    Return the Clock set of CpGs
    
    Parameters
    ----------
    data : dict
        Dictionary returned by getDataDict
    neutral_DNA_CpG_list : list
        Neutral CpGs returned by getNeutralDNACpGs
    n_select : numerical
        How many CpGs to return
    
    Returns
    -------
    Clock_CpGs : ndarray
        Names of CpGs that fit CLOCK_CRITERIA
    """
    
    for cohort in ['tumor', 'basal', 'luminal', 'normals']:
        data[cohort]['beta_values_SELECTION'] = data[cohort]['beta_values'][data[cohort]['pureSamples']]
    
    n_tumors = data['tumor']['beta_values_SELECTION'].shape[1]
    n_basals = data['basal']['beta_values_SELECTION'].shape[1]
    n_luminals = data['luminal']['beta_values_SELECTION'].shape[1]
    n_normals = data['normals']['beta_values_SELECTION'].shape[1]
    print(f"Selecting CpGs with {n_tumors} TCGA tumors, {n_basals} normal basal samples, {n_luminals} normal luminal samples, {n_normals} normal bulk samples")
    
    # Add beta_means, beta_stds, beta_nans to data dict
    addMeanStdsNans(data)
    
    if 'allCpGs' not in data:
        data['allCpGs'] = data['tumor']['beta_values'].index.values
    Clock_CpGs = getCpG_list(data, CLOCK_CRITERIA, starting_CpG_list=neutral_DNA_CpG_list, n_select=n_select)
    
    return Clock_CpGs