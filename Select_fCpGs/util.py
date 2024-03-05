import pandas as pd
import numpy as np
import os
import sys
from MolecularClocks.src.invasiveCpGs_consts import getConsts
from MolecularClocks.src.util import combineFilters
import MolecularClocks.src.methylation_util as m_util

consts = getConsts()
TCGA_datadir = os.path.join(consts['official_indir'], 'TCGA')

# Dictionary defines acceptable ranges (inclusive) for each value
BALANCED_CRITERIA = {
    'normal':{
        'beta_nans':(0, 20),
        'beta_means':(0.4, 0.6)
    },
    'tumor':{
        'beta_nans':(0, 20),
        'beta_means':(0.4, 0.6)
    }
}

def getNeutralDNACpGs():
    TCGA_datadir = os.path.join(consts['official_indir'], 'TCGA')
    
    chip450_info = pd.read_table(os.path.join(TCGA_datadir, 'chip450_annot_cleaned.txt'), index_col=0, low_memory=False)
    chip850_info = pd.read_table(os.path.join(TCGA_datadir, 'chip850_annot_cleaned.txt'), index_col=0, low_memory=False)
    
    neutral_filters = [
        (chip450_info['Regulatory_Feature_Group'].isna()),
        ~chip450_info['HAS_geneName'],
        chip850_info['Regulatory_Feature_Group'].isna(),
        chip850_info['UCSC_RefGene_Name'].isna(),
    ]
    
    neutral_DNA_CpG_list = chip450_info.index[combineFilters(neutral_filters)].values
    
    ## Sanity checks - IMPORTANT - experienced weird errors before with selecting CpGs that did not fit the criteria...
    neutral_chip450_info = chip450_info.loc[neutral_DNA_CpG_list]
    neutral_chip850_info = chip850_info.loc[neutral_DNA_CpG_list]
    assert (neutral_chip450_info['Regulatory_Feature_Group'].isna()).all()
    assert (~neutral_chip450_info['HAS_geneName']).all()
    assert neutral_chip850_info['Regulatory_Feature_Group'].isna().all()
    assert neutral_chip850_info['UCSC_RefGene_Name'].isna().all()

    return neutral_DNA_CpG_list

def getDataDict():
    # Use dictionay to manage two datasets
    # Importing data could take a few minutes once file is downloaded (if on an external file source)
    data = {'tumor':{}, 'normal':{}}
    data['tumor']['beta_values'] = pd.read_table(os.path.join(TCGA_datadir, 'cohort1.methyl.tsv'), index_col=0)
    data['normal']['beta_values'] = consts['project_data_func_dict']['Johnson']['all']()

    # CpG Sites are the same in each DF
    idx_tumor = data['tumor']['beta_values'].index
    idx_normal = data['normal']['beta_values'].index
    assert idx_tumor.isin(idx_normal).all()
    assert idx_normal.isin(idx_tumor).all()

    # Use the same order of CpGs in the index
    data['allCpGs'] = idx_tumor.values
    data['normal']['beta_values'] = data['normal']['beta_values'].loc[data['allCpGs']]

    # Import clinical data
    data['tumor']['clinical'] = pd.read_table(os.path.join(TCGA_datadir, 'cohort1.clinical.tsv'), index_col=0)
    assert data['tumor']['clinical']['submitter_id'].is_unique
    data['tumor']['clinical'] = data['tumor']['clinical'].set_index('submitter_id')

    data['tumor']['ductal_patients'] = data['tumor']['clinical'].index[
        data['tumor']['clinical']['in_analysis_dataset'] & data['tumor']['clinical']['in_CpG_dataset']
    ].values

    # Clip unnecessary accessions from sample IDs
    shortenSampleIDs = lambda x:'-'.join(x.split('-')[:4])
    data['tumor']['beta_values'] = data['tumor']['beta_values'].rename(shortenSampleIDs, axis='columns')
    
    print(f"Selecting CpGs with {data['tumor']['beta_values'].shape[1]} tumors")
    
#     # Import purity estimates
#     # purityEstimates = consts['purityEstimates']()
#     purityEstimates = pd.read_csv(os.path.join(TCGA_datadir, 'ncomms9971-s2.csv'), index_col=0)

#     # Remove duplicate tumors and those without a purity estimate
    
#     duplicated_bool = data['tumor']['beta_values'].columns.duplicated(keep=False)
#     print(f'Removing {duplicated_bool.sum()} tumors for not being from unique patients')
#     print(data['tumor']['beta_values'].columns[duplicated_bool])
#     data['tumor']['beta_values'] = data['tumor']['beta_values'].loc[:, ~duplicated_bool]

#     noPurity_bool = ~data['tumor']['beta_values'].columns.isin(purityEstimates.index)
#     print(f'Removing {noPurity_bool.sum()} tumors for not having a purity estimate')
#     print(data['tumor']['beta_values'].columns[noPurity_bool])
#     data['tumor']['beta_values'] = data['tumor']['beta_values'].loc[:, ~noPurity_bool]
    
#     # Get purity values of remaining tumors/normals
#     data['tumor']['purity'] = purityEstimates.loc[data['tumor']['beta_values'].columns, 'CPE']
#     data['normal']['purity'] = m_util.getLUMP_values(data['normal']['beta_values'])

#     for sample in ['tumor', 'normal']:
#         data[sample]['pureSamples'] = data[sample]['purity'].index[data[sample]['purity'] >= data[sample]['purity_threshold']].values
    
    return data

def addMeanStdsNans(data):
    for sample in ['tumor', 'normal']:
#         if 'beta_values_CpG_SELECTION' not in data[sample]:
# #             data[sample]['beta_values_PURE'] = data[sample]['beta_values'][data[sample]['pureSamples']]

            
#         data[sample]['beta_means'] = data[sample]['beta_values_CpG_SELECTION'].mean(axis=1)
#         data[sample]['beta_stds'] = data[sample]['beta_values_CpG_SELECTION'].std(axis=1)
#         data[sample]['beta_nans'] = data[sample]['beta_values_CpG_SELECTION'].isna().sum(axis=1)

        data[sample]['beta_means'] = data[sample]['beta_values_SELECTION'].mean(axis=1)
        data[sample]['beta_stds'] = data[sample]['beta_values_SELECTION'].std(axis=1)
        data[sample]['beta_nans'] = data[sample]['beta_values_SELECTION'].isna().sum(axis=1)

def gen_CpG_set(data, neutral_DNA_CpG_list, only_ductals=False, n_select=500):
#     for sample in ['tumor', 'normal']:
# #         data[sample]['beta_values_PURE'] = data[sample]['beta_values'][data[sample]['pureSamples']]
    
    data['normal']['beta_values_SELECTION'] = data['normal']['beta_values']
    if only_ductals:
        sampleIDtoPatientID = lambda x:'-'.join(x.split('-')[:-1])
#         ductal_bool = data['tumor']['beta_values_CpG_SELECTION'].columns.to_series().apply(sampleIDtoPatientID).isin(data['tumor']['ductal_patients'])
        ductal_bool = data['tumor']['beta_values'].columns.to_series().apply(sampleIDtoPatientID).isin(data['tumor']['ductal_patients'])
        data['tumor']['beta_values_SELECTION'] = data['tumor']['beta_values'].loc[:, ductal_bool]
    else:
        data['tumor']['beta_values_SELECTION'] = data['tumor']['beta_values']
    
    
    
    # Add beta_means, beta_stds, beta_nans to data dict
    addMeanStdsNans(data)
    
    if 'allCpGs' not in data:
        data['allCpGs'] = data['tumor']['beta_values'].index.values
    balanced_CpGs = m_util.getCpG_list(data, BALANCED_CRITERIA, starting_CpG_list=neutral_DNA_CpG_list, n_select=n_select)
    
    return balanced_CpGs