"""
WCHS_revision.py
=======
Author - Daniel Monyak
5-27-25
=======

Source code for processing the WCHS data

1. Import beta values
2. Filter samples using LUMP purity and eliminate DCIS samples
3. Create table of beta values of Clock sites only for pure samples only

"""

import pandas as pd
import numpy as np
import os
import EpiClockInvasiveBRCA.src.util as epi_util
consts = epi_util.consts


proj_dir = os.path.join(consts['official_indir'], 'WCHS')

# Import calculated beta values
beta_values = pd.read_table(os.path.join(proj_dir, 'beta_values.ClockAndLUMP.manual.transposed.tsv'), index_col=0).T
beta_values = beta_values[np.unique(beta_values.columns)]

# Import sample annotations
sample_annotations = pd.read_csv(os.path.join(proj_dir, 'wchs_450k_pheno_data.csv'), index_col=0)
sample_annotations = sample_annotations.rename(index=lambda x:x[1:]) # eliminate X prefix of sample names
sample_annotations['tumor_stage'] = 'Stage ' + sample_annotations['tumor_stage'].apply(lambda x:None if np.isnan(x) else str(int(x)))

# Exclude DCIS samples
dcis_samples = sample_annotations.index[sample_annotations['tumor_stage'] == 'Stage 0'].values
drop_mask = beta_values.columns.isin(dcis_samples)
print(f'{np.sum(drop_mask)} DCIS (stage 0) samples were dropped')
beta_values = beta_values.loc[:, ~drop_mask]


## Filter samples based on LUMP purity value

LUMP_THRESH = consts['lump_threshold_dict']['Default']

LUMP_purity = epi_util.getLUMP_values(beta_values)
pureSamples = LUMP_purity.index[LUMP_purity >= LUMP_THRESH].to_numpy()
LUMP_purity.to_csv(os.path.join(proj_dir, 'LUMP_purity.txt'), sep='\t')

print(f'{LUMP_purity.shape[0] - len(pureSamples)} samples were dropped for low purity')
print(f'{len(pureSamples)} samples were kept')

## Save beta values of Clock sites for pure samples

Clock_CpGs = np.loadtxt(os.path.join(consts['repo_dir'], '3. Select fCpGs', 'outputs_revision', 'Clock_CpGs.txt'), dtype=str)

beta_values_Clock_CpGs_pureSamples = beta_values.loc[Clock_CpGs, pureSamples]
beta_values_Clock_CpGs_pureSamples.to_csv(os.path.join(proj_dir, 'beta_values_Clock_CpGs_revision_pureSamples.txt'), sep='\t')