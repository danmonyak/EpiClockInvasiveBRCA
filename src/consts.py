"""
consts.py
=======
Author - Daniel Monyak
6-7-24
=======


Provides
    Directories and variables used by other notebooks and source files

"""

import os
import sys

consts = {}

#! Set to parent directory of repo
labdir = os.getenv('lab')
if labdir is None:
    sys.exit('Must set $lab environmental variable to parent directory of MoleculClocks repository...')
    
#! Set to directory with data
consts['official_indir'] = os.path.join(os.getenv("HOME"), 'Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/MolClock_Paper_1/1. Analytic Datasets')
consts['repo_dir'] = os.path.join(labdir, 'EpiClockInvasiveBRCA')
consts['TCGA_datadir'] = os.path.join(consts['official_indir'], 'TCGA')
consts['repo_datadir'] = os.path.join(consts['repo_dir'], 'data')

consts['lump_threshold_dict'] = {
    'Johnson':0.7,
    'Lund':0.6,
    'Reyngold':0.6,
    'Desmedt':0.6
}

consts['palette_jco'] = ["#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#4A6990FF"]

