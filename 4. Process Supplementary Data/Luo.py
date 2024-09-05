"""
Luo.py
=======
Author - Daniel Monyak
9-5-24
=======

Source code for processing the Luo data sample annotations

"""

import pandas as pd
import os
consts = epi_util.consts

# Indir of data
proj_dir = os.path.join(consts['official_indir'], 'Luo')

# Sample Map
sample_map = pd.read_csv(os.path.join(proj_dir, 'sample_map.txt'), sep='\t', names=['GSM', 'sourceName'])
sample_map['Patient'] = sample_map['sourceName'].str.slice(stop=2)
sample_map['Section'] = sample_map['sourceName'].str.slice(start=2)
sample_map = sample_map.drop('sourceName', axis=1)
sample_map.to_csv(os.path.join(proj_dir, 'sample_map_cleaned.txt'), sep='\t', index=False)
