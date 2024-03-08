# import EpiClockInvasiveBRCA.src.cleanAnnotations as cleanAnnot

import pandas as pd
import numpy as np
import os
import sys
from MolecularClocks.src.invasiveCpGs_consts import getConsts
consts = getConsts()

# 450k and 850k microarray site annotations
def chip450_850():
    keep_chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX']
    chip450_info = pd.read_csv(os.path.join(consts['indir'], 'misc', 'cpg_info.csv'), index_col=0)
    float_chrom_rows = chip450_info['CHR'].astype(str).str.contains('\.')
    chip450_info.loc[float_chrom_rows, 'CHR'] = chip450_info.loc[float_chrom_rows, 'CHR'].astype(str).str.split('\.', expand=True)[0]
    chip450_info['CHR'] = 'chr' + chip450_info['CHR'].astype(str)
    chip450_info = chip450_info.loc[chip450_info['CHR'].isin(keep_chroms)]
    chip450_info = chip450_info.loc[chip450_info.index.str.startswith('cg')]
    assert (chip450_info['MAPINFO'] == chip450_info['MAPINFO'].astype(int)).all()
    chip450_info['MAPINFO'] = chip450_info['MAPINFO'].astype(int)
    chip450_info['Genome_Build'] = chip450_info['Genome_Build'].astype(int).astype(str)
    
    
    assert (chip450_info['Relation_to_UCSC_CpG_Island'] == 'None').sum() == 0
    assert (chip450_info['Regulatory_Feature_Group'] == 'None').sum() == 0
    chip450_info['Relation_to_UCSC_CpG_Island'] = chip450_info['Relation_to_UCSC_CpG_Island'].fillna('None')
    chip450_info['Regulatory_Feature_Group'] = chip450_info['Regulatory_Feature_Group'].fillna('None')
    chip450_info['HAS_geneName'] = ~(chip450_info['UCSC_RefGene_Name'].isna())
    #####
    #####
    assert not chip450_info.index.duplicated().any()

    chip850_info = pd.read_csv(os.path.join(consts['indir'], 'misc', 'chip850_info.csv'), index_col='Name')
    chip850_info = chip850_info.set_index('Methyl450_Loci')
    chip850_info = chip850_info.loc[chip850_info.index.isin(chip450_info.index)]
    chip850_info = chip850_info.loc[~(chip850_info.index.duplicated())]
    chip450_info = chip450_info.loc[chip850_info.index]

    chip450_info.to_csv(os.path.join(consts['official_indir'], 'TCGA', 'chip450_annot_cleaned.txt'), sep='\t')
    chip850_info.to_csv(os.path.join(consts['official_indir'], 'TCGA', 'chip850_annot_cleaned.txt'), sep='\t')

def TCGA_clinical():
    TCGA_datadir = os.path.join(consts['official_indir'], 'TCGA')
    
    clinical = pd.read_table(os.path.join(TCGA_datadir, 'cohort1.clinical.tsv'), index_col=0, na_values="'--")
    
#     clinical = clinical.loc[~clinical['ajcc_pathologic_t'].isna()]
    
    ajcc_pathologic_t = clinical['ajcc_pathologic_t'].apply(lambda x:x[:2])
    ajcc_pathologic_stage = clinical['ajcc_pathologic_stage'].str.rstrip('ABC')
    
    clinical = clinical.drop(['ajcc_pathologic_t', 'ajcc_pathologic_stage'], axis=1)
    clinical['ajcc_pathologic_t'] = ajcc_pathologic_t
    clinical['ajcc_pathologic_stage'] = ajcc_pathologic_stage
    return clinical
#     clinical.to_csv(os.path.join(TCGA_datadir, 'cohort1.clinical.cleaned.tsv'), sep='\t')

