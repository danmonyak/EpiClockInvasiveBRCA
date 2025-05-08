"""
Aurora.py
=======
Author - Daniel Monyak
9-5-24
=======

Source code for processing data from the AURORA US Metastasis project

1. Process sample annotations
2. Create table of beta values of Clock sites only
3. Filter samples using LUMP purity
4. Pick primary and metastasis for each patient
5. Make sample map

"""

import pandas as pd
import numpy as np
import os
import EpiClockInvasiveBRCA.src.util as epi_util
consts = epi_util.consts

proj_dir = os.path.join(consts['official_indir'], 'Aurora')

# Sample_annotations

sample_annotations = pd.read_excel(os.path.join(proj_dir, '43018_2022_491_MOESM2_ESM', 'Supplementary_Table.2.xlsx'), sheet_name='2.AURORA study')
sample_annotations.index = sample_annotations['BCR Portion barcode'].map(lambda x:'-'.join(x.split('.')))

## Beta values

beta_values = pd.read_table(os.path.join(proj_dir, 'GSE212370_20220828_MethylationEPIC_GEO_data_matrix146.tsv'), index_col=0)
beta_values = beta_values.rename(columns=lambda x:'-'.join(x.split('-')[:6]))

Clock_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], '3. Select fCpGs', 'outputs_revision', 'Clock_CpGs.txt'), dtype=str)
# Some might be missing
before = Clock_CpGs.shape[0]
Clock_CpGs = np.intersect1d(Clock_CpGs, beta_values.index)
after = Clock_CpGs.shape[0]
print(f'{before - after} Clock CpGs were not measured in the Aurora dataset')
#
beta_values_Clock = beta_values.loc[Clock_CpGs]

## LUMP purity - sample filtering
## Filter samples based on LUMP purity value

LUMP_THRESH = consts['lump_threshold_dict']['Default']

LUMP_purity = epi_util.getLUMP_values(beta_values)
sample_annotations['LUMP'] = LUMP_purity
sample_annotations['pure'] = sample_annotations['LUMP'] >= LUMP_THRESH
LUMP_purity.to_csv(os.path.join(proj_dir, 'LUMP_purity.txt'), sep='\t')

## Pick one primary and one metastasis for each patient
## Only pick patients with one pure sample in each category
## Pick highest LUMP sample if there are more than one

patient_list = []
prim_list = []
met_list = []
patient_multiple_mets = {}

all_patients = sample_annotations['Patient'].unique()
print(f'Initially using {all_patients.shape[0]} patients')
for patient in all_patients:
    
    sa_patient = sample_annotations.loc[(sample_annotations['Patient'] == patient) & sample_annotations['pure']]
    
    prim_samps = sa_patient.loc[sa_patient['Sample Type'] == 'Primary', 'LUMP'].sort_values(ascending=False).index.values
    met_samps = sa_patient.loc[sa_patient['Sample Type'] == 'Metastasis', 'LUMP'].sort_values(ascending=False).index.values
    
    # If more than one to choose from
    if (prim_samps.shape[0] > 0) and (met_samps.shape[0] > 0):
#         met_samps = np.abs(sa_patient.loc[met_samps, 'LUMP'] - sa_patient.loc[prim_samps[0], 'LUMP']).sort_values(ascending=True).index.values
        patient_list.append(patient)
        prim_list.append(prim_samps[0])
        met_list.append(met_samps[0])

    if met_samps.shape[0] > 1:
        patient_multiple_mets[patient] = met_samps.tolist()

print(f'Final cohort: {len(patient_list)} patients')

# Compile into sample_map

patient_map = pd.DataFrame({'Patient':patient_list, 'Primary':prim_list, 'Metastasis':met_list})
sample_map = patient_map.melt(id_vars='Patient').rename(columns={'variable':'Sample Type'}).set_index('value')
sample_map = sample_map.sort_values(by=['Patient', 'Sample Type'], ascending=[True, False])
sample_map.to_csv(os.path.join(proj_dir, 'sample_map.txt'), sep='\t')

## Save beta values of Clock sites for pure samples

beta_values_Clock_CpGs_pureSamples = beta_values_Clock[sample_annotations.index[sample_annotations['pure']]]
beta_values_Clock_CpGs_pureSamples.to_csv(os.path.join(proj_dir, 'beta_values_Clock_CpGs_revision_pureSamples.txt'), sep='\t')

## Save list of patients that has a pure primary and metastasis

np.savetxt(os.path.join(proj_dir, 'patients_bothSamplesPure.txt'), patient_list, fmt='%s')

## Save list of metastasis per patient with multiple metastases
pd.Series(patient_multiple_mets).to_csv(os.path.join(proj_dir, 'patient_multiple_mets.txt'), sep='\t')