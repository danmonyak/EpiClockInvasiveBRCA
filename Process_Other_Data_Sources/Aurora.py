import pandas as pd
import numpy as np
import os
import EpiClockInvasiveBRCA.src.util as epi_util
from EpiClockInvasiveBRCA.src.consts import consts

proj_dir = os.path.join(consts['official_indir'], 'Aurora')

# sample map

sample_annotations = pd.read_excel(os.path.join(proj_dir, '43018_2022_491_MOESM2_ESM', 'Supplementary_Table.2.xlsx'), sheet_name='2.AURORA study')
sample_annotations.index = sample_annotations['BCR Portion barcode'].map(lambda x:'-'.join(x.split('.')))

## Beta values

beta_values = pd.read_table(os.path.join(proj_dir, 'GSE212370_20220828_MethylationEPIC_GEO_data_matrix146.tsv'), index_col=0)
beta_values = beta_values.rename(columns=lambda x:'-'.join(x.split('-')[:6]))

Clock_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], 'Select_fCpGs', 'outputs', 'balanced_CpGs.txt'), dtype=str)
beta_values_Clock = beta_values.loc[Clock_CpGs]

## LUMP purity - sample filtering

LUMP_THRESH = 0.6

LUMP_purity = epi_util.getLUMP_values(beta_values)
sample_annotations['LUMP'] = LUMP_purity
sample_annotations['pure'] = sample_annotations['LUMP'] >= LUMP_THRESH
LUMP_purity.to_csv(os.path.join(proj_dir, 'LUMP_purity.txt'), sep='\t')

## Pick primary and metastasis for each patient

patient_list = []
prim_list = []
met_list = []

for patient in sample_annotations['Patient'].unique():
    sa_patient = sample_annotations.loc[(sample_annotations['Patient'] == patient) & sample_annotations['pure']]
    
    prim_samps = sa_patient.loc[sa_patient['Sample Type'] == 'Primary', 'LUMP'].sort_values(ascending=False).index.values
    met_samps = sa_patient.loc[sa_patient['Sample Type'] == 'Metastasis', 'LUMP'].sort_values(ascending=False).index.values
    
    if (prim_samps.shape[0] > 0) and (met_samps.shape[0] > 0):
#         met_samps = np.abs(sa_patient.loc[met_samps, 'LUMP'] - sa_patient.loc[prim_samps[0], 'LUMP']).sort_values(ascending=True).index.values
        patient_list.append(patient)
        prim_list.append(prim_samps[0])
        met_list.append(met_samps[0])

patient_map = pd.DataFrame({'Patient':patient_list, 'Primary':prim_list, 'Metastasis':met_list})
sample_map = patient_map.melt(id_vars='Patient').rename(columns={'variable':'Sample Type'}).set_index('value')
sample_map = sample_map.sort_values(by=['Patient', 'Sample Type'], ascending=[True, False])
sample_map.to_csv(os.path.join(proj_dir, 'sample_map.txt'), sep='\t')

## Save beta values of Clock sites for pure samples

beta_values_Clock_CpGs_pureSamples = beta_values_Clock[sample_annotations.index[sample_annotations['pure']]]
beta_values_Clock_CpGs_pureSamples.to_csv(os.path.join(proj_dir, 'beta_values_Clock_CpGs_pureSamples.txt'), sep='\t')

## Save list of patients that has a pure primary and metastasis

np.savetxt(os.path.join(proj_dir, 'patients_bothSamplesPure.txt'), patient_list, fmt='%s')