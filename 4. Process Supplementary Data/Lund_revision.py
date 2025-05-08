"""
Lund.py
=======
Author - Daniel Monyak
9-5-24
=======

Source code for processing the Lund cohort data

1. Process sample annotations
2. Import beta values
3. Filter samples
4. Calculate c_beta values for samples

"""

import pandas as pd
import numpy as np
import os
import EpiClockInvasiveBRCA.src.util as epi_util
consts = epi_util.consts

proj_dir = os.path.join(consts['official_indir'], 'Lund')

## Create outdir if necessary
output_dir = os.path.join('outputs_revision', 'Lund')
os.makedirs(output_dir, exist_ok=True)

## Sample annotations

# File mapping TAX IDs to GSM IDs
TAX_to_GSM_mapping = pd.read_table(os.path.join(proj_dir, 'TAX_to_GSM_mapping.txt'), index_col=0).squeeze('columns')

# Clinical file
GSE75067_clinical = pd.read_table(os.path.join(proj_dir, 'GSE75067_sample_annotations.txt'), index_col=0)
clinical = GSE75067_clinical.reset_index().merge(TAX_to_GSM_mapping, left_on='Title', right_index=True).set_index('GSM')
clinical['grade'] = clinical['grade'].map({1:'Grade 1', 2:'Grade 2', 3:'Grade 3'}).fillna(value='None')
###################

print('Loading beta values, should take <3 minutes...')

## Beta values
beta_values = pd.read_table(os.path.join(proj_dir, 'GSE75067_betaValues.txt'), index_col=0)

print('Loaded.')

print(f'{clinical.shape[0]} tumors initially')

# Remove second tumor from patient with two tumors
# Chose the tumor that would have been excluded anyway
# 1 was a second primary with mixed histology, 2 were metastases, 1 was lobular histology
clinical = clinical.drop(['GSM1941866', 'GSM1941946', 'GSM1942009', 'GSM1941877'], axis=0)
print(f'Manualy dropped 4 tumors from patient with multiple tumors')

# Only include primary tumors
clinical['in_CpG_dataset'] = clinical['primary'] == 'primary'
print(f'{clinical["in_CpG_dataset"].sum()} primary tumors')
print(f'{clinical.loc[clinical["in_CpG_dataset"], "age"].unique().shape[0]} unique patients')

## Filter samples based on LUMP purity value
# Set in_CpG_dataset and reason_purity flags accordingly - based on LUMP value
# Remove impure tumors from beta_values

# LUMP threshold for this dataset
LUMP_THRESH = consts['lump_threshold_dict']['Lund']

LUMP_purity = epi_util.getLUMP_values(beta_values)
clinical = clinical.merge(LUMP_purity.rename('LUMP'), left_index=True, right_index=True)
clinical['reason_purity'] = clinical['in_CpG_dataset'] & (clinical['LUMP'] < LUMP_THRESH)
clinical['in_CpG_dataset'] &= ~clinical['reason_purity']
print(f'{clinical["reason_purity"].sum()} tumors removed for having LUMP < {LUMP_THRESH}')

# revision sites
Clock_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], '3. Select fCpGs', 'outputs_revision', 'Clock_CpGs.txt'), dtype=str)

# Don't need to remove any samples for having >= 5% missing values in Clock fCpGs
assert (beta_values.loc[Clock_CpGs].isna().mean(axis=0) < 0.05).all()
print('0 tumors had to be excluded for too many missing Clock beta values')

beta_values = beta_values[clinical.index[clinical['in_CpG_dataset']]]
print(f'{clinical["in_CpG_dataset"].sum()} tumors kept')

# Set in_analysis_dataset flag to indicate that tumor is ductal or not ductal
clinical['in_analysis_dataset'] = clinical['in_CpG_dataset'] & (clinical['tumorType_loRes']=='ductal')
n_removed_histology = clinical['in_CpG_dataset'].sum() - clinical['in_analysis_dataset'].sum()
print(f'{n_removed_histology} tumors removed for having non-ductal histology')
print(f'{clinical["in_analysis_dataset"].sum()} final tumors')

# Save clinical file
clinical.to_csv(os.path.join(proj_dir, 'Lund.clinical.revision.txt'), sep='\t')

print('Saved modified clinical file.')

## Calculate and save c_beta values for each tumor
beta_values_Clock_CpGs_pureSamples = beta_values.loc[Clock_CpGs, clinical.index[clinical['in_analysis_dataset']]]
c_beta = 1 - beta_values_Clock_CpGs_pureSamples.std(axis=0)
c_beta.name = 'c_beta'
c_beta.to_csv(os.path.join(output_dir, 'Lund.c_beta.txt'), sep='\t')

print('Saved c_beta value of each tumor.')


## Save beta values of Clock sites for pure samples
beta_values_Clock_CpGs_pureSamples.to_csv(os.path.join(proj_dir, 'beta_values_Clock_CpGs_revision_pureSamples.tsv'), sep='\t')

