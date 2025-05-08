"""
Germany_revision.py
=======
Author - Daniel Monyak
4-11-25
=======

Source code for processing the Germany data

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


proj_dir = os.path.join(consts['official_indir'], 'Germany')

## Sample annotations

sample_annotations = pd.read_table(os.path.join(proj_dir, 'sample_annotations.txt'), index_col=0, na_values=['NA'])

# Process all columns
data_dict = {}
for i in range(sample_annotations.shape[1]):
    col = sample_annotations.iloc[:, i]
    for j in range(col.shape[0]):
        if (col.index[j] != 'Sample_characteristics_ch1') or (type(col.iloc[j]) is float):
            continue
        label, data = col.iloc[j].split(': ')
        if label not in data_dict:
            data_dict[label] = {}
        data_dict[label][col.name] = data

sample_annotations_clean = pd.DataFrame(data_dict)

# Add clean status column
status_mapper = {0:'normal', 1:'normal-adjacent', 2:'breast cancer', 3:'normal-brca1', 4:'cancer-brca1'}
sample_annotations_clean['status'] = sample_annotations_clean['status(0=normal 1=normal-adjacent 2=breast cancer 3=normal-brca1 4=cancer-brca1)'].astype(int).map(status_mapper)

####################################

## Beta values
beta_values = pd.read_table(os.path.join(proj_dir, 'beta_values.txt'), index_col=0)

## Filter samples based on LUMP purity value

LUMP_THRESH = consts['lump_threshold_dict']['Default']

LUMP_purity = epi_util.getLUMP_values(beta_values)
pureSamples = LUMP_purity.index[LUMP_purity >= LUMP_THRESH].to_numpy()
LUMP_purity.to_csv(os.path.join(proj_dir, 'LUMP_purity.txt'), sep='\t')

## Save beta values of Clock sites for pure samples

Clock_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], '3. Select fCpGs', 'outputs_revision', 'Clock_CpGs.txt'), dtype=str)
beta_values_Clock_CpGs_pureSamples = beta_values.loc[Clock_CpGs, pureSamples]
beta_values_Clock_CpGs_pureSamples.to_csv(os.path.join(proj_dir, 'beta_values_Clock_CpGs_revision_pureSamples.txt'), sep='\t')

## Process other columns of sample_annotations_clean
sample_annotations_clean['Ki-67'] = sample_annotations_clean['ki67'].map({'0':'Ki-67 Negative', '1':'Ki-67 Positive'})
sample_annotations_clean['c_beta'] = 1 - beta_values_Clock_CpGs_pureSamples.std(axis=0, numeric_only=True)
sample_annotations_clean['LUMP'] = LUMP_purity

## Restrict to breast cancer cohort and pure samples
sample_annotations_clean['in_analysis_dataset'] = (sample_annotations_clean['status'] == 'breast cancer') & sample_annotations_clean.index.isin(beta_values_Clock_CpGs_pureSamples.columns)

## Save file
sample_annotations_clean.to_csv(os.path.join(proj_dir, 'sample_annotations_clean.txt'), sep='\t')
