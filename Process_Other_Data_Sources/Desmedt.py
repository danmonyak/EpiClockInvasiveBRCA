import pandas as pd
import numpy as np
import os
import EpiClockInvasiveBRCA.src.util as epi_util
from MolecularClocks.src.invasiveCpGs_consts import getConsts
consts = getConsts()

proj_dir = os.path.join(consts['official_indir'], 'Desmedt')

# sample map

sample_annotations = pd.read_table(os.path.join(proj_dir, 'sample_annotations.txt'), index_col=0).T

col_list = []
for i in range(sample_annotations.shape[1]):
    col = sample_annotations.iloc[:, i]
    if col.name == 'Sample_characteristics_ch1':
        col_split = col.str.split(': ', expand=True)
        assert col_split[0].unique().shape[0] == 1
        col_list.append(col_split[1].rename(col_split.iloc[0, 0]))
    else:
        col_list.append(col)
sample_annotations_clean = pd.concat(col_list, axis=1)

sample_annotations_clean['foci'] = sample_annotations_clean['Sample_description'].apply(lambda x:x[-1]).astype(str)
ser = sample_annotations_clean['Sample_description']
ser = ser.str.rstrip('_12')
sample_annotations_clean['patient'] = ser.str.split('_foci', expand=True)[0].apply(lambda x:x[-1])
sample_annotations_clean = sample_annotations_clean.drop('Sample_description', axis=1)

sample_counts = sample_annotations_clean.groupby('patient')['foci'].value_counts().unstack()
assert (sample_counts == 1).all(axis=None)

sample_annotations_clean.to_csv(os.path.join(proj_dir, 'sample_annotations_clean.txt'), sep='\t')

####################################

# pure samples

LUMP_THRESH = 0.6

beta_values = pd.read_table(os.path.join(proj_dir, 'beta_values.txt'), index_col=0)

LUMP_purity = epi_util.getLUMP_values(beta_values)
pureSamples = LUMP_purity.index[LUMP_purity >= LUMP_THRESH].to_numpy()
# np.savetxt(os.path.join(proj_dir, 'pureSamples.txt'), pureSamples, fmt='%s')

balanced_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], 'Select_fCpGs', 'outputs', 'balanced_CpGs.txt'), dtype=str)

beta_values_balanced_CpGs_pureSamples = beta_values.loc[balanced_CpGs, pureSamples]
beta_values_balanced_CpGs_pureSamples.to_csv(os.path.join(proj_dir, 'beta_values_balanced_CpGs_pureSamples.txt'), sep='\t')

# sample_map_with_lump = pd.concat([sample_annotations_clean, LUMP_purity.rename('lump')], axis=1)
# pureSample_count = sample_map_with_lump.loc[pureSamples, 'patient'].value_counts()
pureSample_count = sample_annotations_clean.loc[pureSamples, 'patient'].value_counts()
patients_bothSamplesPure = pureSample_count.index[pureSample_count == 2].values
np.savetxt(os.path.join(proj_dir, 'patients_bothSamplesPure.txt'), patients_bothSamplesPure, fmt='%s')