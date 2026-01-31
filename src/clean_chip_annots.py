"""
clean_chip_annots.py
=======
Author - Daniel Monyak
10-1-24
=======

Script that needs to be run to process the 450k and 850k array annotation files
"""

import pandas as pd
import numpy as np
import os
import sys
import EpiClockInvasiveBRCA.src.util as epi_util
consts = epi_util.consts

# Directory with annotation files
chip_annot_dir = os.path.join(consts['TCGA_datadir'], 'chip_annots')

# Chromosome lists
keep_chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

## Load raw annotation files
chip450_info = pd.read_csv(os.path.join(chip_annot_dir, 'humanmethylation450_15017482_v1-2.csv'), index_col=0, low_memory=False)
print('Loaded 450K annotation')

chip850_info = pd.read_csv(os.path.join(chip_annot_dir, 'infinium-methylationepic-v-1-0-b5-manifest-file.csv'), index_col=0, low_memory=False)
print('Loaded 850K annotation')

###

## Process 450k annotations
chip450_info['CHR'] = 'chr' + chip450_info['CHR'].astype(str)                  # convert chromosome column to str
chip450_info = chip450_info.loc[chip450_info['CHR'].isin(keep_chroms)]         # remove CpGs not on chromosomes
chip450_info = chip450_info.loc[chip450_info.index.str.startswith('cg')]       # remove markers that are not CpGs
assert (chip450_info['MAPINFO'] == chip450_info['MAPINFO'].astype(int)).all()  # sanity check that no info will be lost on conversion
chip450_info['MAPINFO'] = chip450_info['MAPINFO'].astype(int)                  # convert coordinate column to int

###

# Sanity check -- no duplicated CpGs
assert not chip450_info.index.duplicated().any()

# Process 850k annotations dataframe
# Align to 450k annotations dataframe
chip850_info = chip850_info.loc[chip850_info['Methyl450_Loci'].fillna(False)]   # Select 850k only sites that are on the 450k
chip850_info = chip850_info.loc[chip850_info.index.isin(chip450_info.index)]    # Select 850k only sites that are on the 450k - check 2
chip850_info = chip850_info.loc[~(chip850_info.index.duplicated())]             # Remove duplicated rows
chip450_info = chip450_info.loc[chip850_info.index]                             # Align 450k dataframe -- same order as 850k dataframe

# Save annotation dataframes
chip450_info.to_csv(os.path.join(chip_annot_dir, 'chip450_annot_cleaned.txt'), sep='\t')
chip850_info.to_csv(os.path.join(chip_annot_dir, 'chip850_annot_cleaned.txt'), sep='\t')
