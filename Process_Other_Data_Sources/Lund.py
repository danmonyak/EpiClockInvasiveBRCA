import pandas as pd
import numpy as np
import os
import MolecularClocks.src.methylation_util as m_util
from MolecularClocks.src.invasiveCpGs_consts import getConsts

consts = getConsts()

# LUMP threshold for this dataset
LUMP_THRESH = 0.6

# Indir of data
proj_dir = os.path.join(consts['official_indir'], 'Ringner')

# File mapping TAX IDs to GSM IDs
TAX_to_GSM_mapping = pd.read_table(os.path.join(proj_dir, 'TAX_to_GSM_mapping.txt'), index_col=0).squeeze('columns')

# Clinical file
GSE75067_clinical = pd.read_table(os.path.join(proj_dir, 'GSE75067_sample_annotations.txt'), index_col=0)
clinical = GSE75067_clinical.reset_index().merge(TAX_to_GSM_mapping, left_on='Title', right_index=True).set_index('GSM')
clinical['grade'] = clinical['grade'].map({1:'Grade 1', 2:'Grade 2', 3:'Grade 3'}).fillna(value='None')
###################

print('Loading beta values, should take <3 minutes...')

# Import beta values
beta_values = pd.read_table(os.path.join(proj_dir, 'GSE75067_betaValues.txt'), index_col=0)

print('Loaded.')

# Create outdir if necessary
output_dir = os.path.join('outputs', 'cohort_T2')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
print(f'{clinical.shape[0]} tumors initially')

# Filter by LUMP
# Set in_CpG_dataset and reason_purity flags accordingly - based on LUMP value
# Remove impure tumors from beta_values
LUMP_purity = m_util.getLUMP_values(beta_values)
clinical = clinical.merge(LUMP_purity.rename('LUMP'), left_index=True, right_index=True)
clinical['reason_purity'] = clinical['LUMP'] < LUMP_THRESH
clinical['in_CpG_dataset'] = ~clinical['reason_purity']
print(f'{clinical["reason_purity"].sum()} tumors removed for having LUMP < {LUMP_THRESH}')

clinical['reason_primary'] = clinical['primary'] != 'primary'
clinical['in_CpG_dataset'] &= ~clinical['reason_primary']
print(f'{clinical["reason_primary"].sum()} tumors removed for not being primary tumors')

filt = clinical.loc[clinical['in_CpG_dataset']]
clinical = clinical.merge(filt['age'].duplicated().rename('reason_duplicate'), left_index=True, right_index=True)
clinical['reason_duplicate'] = clinical['reason_duplicate'].fillna(False)
clinical['in_CpG_dataset'] &= ~clinical['reason_duplicate']
print(f'{clinical["reason_duplicate"].sum()} tumors removed for being from the same patient as a previous one.')

beta_values = beta_values[clinical.index[clinical['in_CpG_dataset']]]
print(f'{clinical["in_CpG_dataset"].sum()} tumors kept')

# Set in_analysis_dataset flag to indicate that tumor is ductal or not ductal
clinical['in_analysis_dataset'] = clinical['in_CpG_dataset'] & (clinical['tumorType_loRes']=='ductal')
clinical.to_csv(os.path.join(proj_dir, 'cohort.T2.clinical.txt'), sep='\t')

print(f'{clinical["in_analysis_dataset"].sum()} ductal tumors kept for analysis')
print('Saved modified clinical file.')

# Calculate and save c_beta values for each tumor
balanced_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], 'Select_fCpGs', 'outputs', 'balanced_CpGs.txt'), dtype=str)
c_beta = 1 - beta_values.loc[balanced_CpGs, clinical.index[clinical['in_analysis_dataset']]].std(axis=0)
c_beta.name = 'c_beta'
c_beta.to_csv(os.path.join(output_dir, 'cohort.T2.c_beta.txt'), sep='\t')

print('Saved c_beta value of each tumor.')

# Don't need to remove any samples for having >= 5% missing values in Clock fCpGs
assert (beta_values.loc[balanced_CpGs].isna().mean(axis=0) < 0.05).all()