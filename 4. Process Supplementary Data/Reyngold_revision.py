"""
Reyngold.py
=======
Author - Daniel Monyak
9-5-24
=======

Source code for processing the Reyngold data

1. Process sample annotations
2. Import beta values
3. Filter samples using LUMP purity
4. Create table of beta values of Clock sites only for pure samples only
5. Pick primary and metastasis for each patient

"""

import pandas as pd
import numpy as np
import os
import EpiClockInvasiveBRCA.src.util as epi_util
consts = epi_util.consts

proj_dir = os.path.join(consts['official_indir'], 'Reyngold')

## Sample annotations

sample_annotations = pd.read_table(os.path.join(proj_dir, 'sample_annotations.txt'), index_col=0).T

# Process all columns
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

sample_counts = sample_annotations_clean.groupby('unique patient id')['sample type'].value_counts().unstack()
assert (sample_counts == 1).all(axis=None)   # Should only be one primary and one metastasis per patient
assert (sample_annotations_clean['tissue'] == 'breast tumor').all()

sample_map = sample_annotations_clean[['unique patient id', 'sample type']]
sample_map.to_csv(os.path.join(proj_dir, 'sample_map.txt'), sep='\t')

####################################

## Beta values

beta_values = pd.read_table(os.path.join(proj_dir, 'GSE58999_betaValues.txt'), index_col=0)

## Filter samples based on LUMP purity value

LUMP_THRESH = consts['lump_threshold_dict']['Default']

LUMP_purity = epi_util.getLUMP_values(beta_values)
pureSamples = LUMP_purity.index[LUMP_purity >= LUMP_THRESH].to_numpy()
LUMP_purity.to_csv(os.path.join(proj_dir, 'LUMP_purity.txt'), sep='\t')

## Save beta values of Clock sites for pure samples

Clock_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], '3. Select fCpGs', 'outputs_revision', 'Clock_CpGs.txt'), dtype=str)
beta_values_Clock_CpGs_pureSamples = beta_values.loc[Clock_CpGs, pureSamples]
beta_values_Clock_CpGs_pureSamples.to_csv(os.path.join(proj_dir, 'beta_values_Clock_CpGs_revision_pureSamples.txt'), sep='\t')

## Save list of patients wth both samples pure

pureSample_count = sample_map.loc[pureSamples, 'unique patient id'].value_counts()
patients_bothSamplesPure = pureSample_count.index[pureSample_count == 2].values
np.savetxt(os.path.join(proj_dir, 'patients_bothSamplesPure.txt'), patients_bothSamplesPure, fmt='%s')