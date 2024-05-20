import pandas as pd
import os
from EpiClockInvasiveBRCA.src.consts import consts

# Indir of data
proj_dir = os.path.join(consts['official_indir'], 'Luo')

# Sample Map
sample_map = pd.read_csv(os.path.join(proj_dir, 'sample_map.txt'), sep='\t', names=['GSM', 'sourceName'])
sample_map['Patient'] = sample_map['sourceName'].str.slice(stop=2)
sample_map['Section'] = sample_map['sourceName'].str.slice(start=2)
sample_map = sample_map.drop('sourceName', axis=1)
sample_map.to_csv(os.path.join(proj_dir, 'sample_map_cleaned.txt'), sep='\t', index=False)
