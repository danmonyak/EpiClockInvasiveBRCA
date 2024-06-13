import pandas as pd
import numpy as np
import os
import EpiClockInvasiveBRCA.src.util as epi_util
from EpiClockInvasiveBRCA.src.consts import consts

TCGA_clinical_dir = os.path.join(consts['official_indir'], 'TCGA')
patient_to_sample_IDs = pd.read_table(os.path.join(TCGA_clinical_dir, 'patient_to_sample_IDs.txt'), index_col=0).squeeze('columns')
sample_IDs_used = patient_to_sample_IDs.values

filepath = '/Users/danielmonyak/Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/MolClock_Paper_1/Data/pan-cancer_TCGA/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv'
rowname_list = []
data_list = []

with open(filepath, 'r') as f:
	header = f.readline()
	specific_sample_list = header.rstrip('\n').replace('"', '').split('\t')[1:]
	sample_list = ['-'.join(s.split('-')[:4]) for s in specific_sample_list]
	
	use_sample = np.isin(sample_list, sample_IDs_used)
	lower = np.min(np.where(use_sample))
	upper = np.max(np.where(use_sample))

	i = 0
	while True:
		if i % 1000 == 0:
			print(i)
		# if i == 100:
		# 	break

		i += 1

		line = f.readline()
		if line == '':
			break

		row_list = line.rstrip('\n').replace('"', '').split('\t')
		rowname = row_list[0]
		data = np.array(row_list[1:])[use_sample]

		rowname_list.append(rowname)
		data_list.append(data)


data_selected = pd.DataFrame(columns=sample_IDs_used, index=rowname_list, data=np.stack(data_list, axis=0))
data_selected.to_csv(os.path.join(TCGA_clinical_dir, 'gene_expr.txt'), sep='\t')
