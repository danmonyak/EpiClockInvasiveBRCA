import pandas as pd
import numpy as np
import os
import sys
import EpiClockInvasiveBRCA.src.util as epi_util
consts = epi_util.consts

## 450k and 850k microarray site annotations
# def chip450_850():
    
chip_annot_dir = os.path.join(consts['TCGA_datadir'], 'chip_annots')
keep_chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX']

chip450_info = pd.read_csv(os.path.join(chip_annot_dir, 'humanmethylation450_15017482_v1-2.csv'), index_col=0, low_memory=False)
print('Loaded 450K annotation')
#     chip850_info = pd.read_csv(os.path.join(consts['indir'], 'misc', 'chip850_info.csv'), index_col='Name', low_memory=False)
#     chip850_info = chip850_info.set_index('Methyl450_Loci')
chip850_info = pd.read_csv(os.path.join(chip_annot_dir, 'infinium-methylationepic-v-1-0-b5-manifest-file.csv'), index_col=0, low_memory=False)
print('Loaded 850K annotation')



## when we were using the old chip450_info file - not necessary with new one
# float_chrom_rows = chip450_info['CHR'].astype(str).str.contains('\.')
# chip450_info.loc[float_chrom_rows, 'CHR'] = chip450_info.loc[float_chrom_rows, 'CHR'].astype(str).str.split('\.', expand=True)[0]
chip450_info['CHR'] = 'chr' + chip450_info['CHR'].astype(str)
chip450_info = chip450_info.loc[chip450_info['CHR'].isin(keep_chroms)]
chip450_info = chip450_info.loc[chip450_info.index.str.startswith('cg')]
assert (chip450_info['MAPINFO'] == chip450_info['MAPINFO'].astype(int)).all()
chip450_info['MAPINFO'] = chip450_info['MAPINFO'].astype(int)
# chip450_info['Genome_Build'] = chip450_info['Genome_Build'].astype(int).astype(str)


# assert (chip450_info['Relation_to_UCSC_CpG_Island'] == 'None').sum() == 0
# assert (chip450_info['Regulatory_Feature_Group'] == 'None').sum() == 0
# chip450_info['Relation_to_UCSC_CpG_Island'] = chip450_info['Relation_to_UCSC_CpG_Island'].fillna('None')
# chip450_info['Regulatory_Feature_Group'] = chip450_info['Regulatory_Feature_Group'].fillna('None')

#####
#####
assert not chip450_info.index.duplicated().any()


chip850_info = chip850_info.loc[chip850_info['Methyl450_Loci'].fillna(False)]
chip850_info = chip850_info.loc[chip850_info.index.isin(chip450_info.index)]
chip850_info = chip850_info.loc[~(chip850_info.index.duplicated())]
chip450_info = chip450_info.loc[chip850_info.index]

chip450_info.to_csv(os.path.join(chip_annot_dir, 'chip450_annot_cleaned.txt'), sep='\t')
chip850_info.to_csv(os.path.join(chip_annot_dir, 'chip850_annot_cleaned.txt'), sep='\t')

print('Done!')