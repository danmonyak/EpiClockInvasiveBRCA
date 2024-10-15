"""
run_ensemble.py
=======
Author - Daniel Monyak
9-4-24
=======

Source code for running a few different variations of the fCpG Ensemble simulation
3 different parameter sets are available
Paramter set is selected in the command line
    python run_ensemble.py [0/1/2]
    
"""

import numpy as np
import sys
import os
from time import process_time
import simulation as sim

# No index provided
if len(sys.argv) == 1:
    sys.exit('Enter index of parameter dict...')

params_idx = int(sys.argv[1])
params_list = [
    {'output_dir':'3_sites', 'n_CpGs_each':1, 'death_rate':0.15, 'nyears':2, 'seed':0},
    {'output_dir':'90_sites', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':0},
    {'output_dir':'larger_death_rate', 'n_CpGs_each':1, 'death_rate':0.16, 'nyears':4, 'seed':0},
]
prog_params = params_list[params_idx]

# Constant parameters
FLIP_RATE = 0.002
PROLIF_RATE = 0.17

os.makedirs(prog_params['output_dir'], exist_ok=True)

init_params = {'flip_rate':FLIP_RATE, # flip rate per cell division per allele
               'prolif_rate':PROLIF_RATE, # cell divisions per day
               'death_rate':prog_params['death_rate'], # cell deaths per day,
               'init_site_state_counts':[prog_params['n_CpGs_each'], prog_params['n_CpGs_each'], 0, prog_params['n_CpGs_each']],
              }
gen = np.random.default_rng(prog_params['seed'])
ensmbl = sim.Ensemble(init_params, gen)

total_days = int(prog_params['nyears'] * 365)

beta_list = []   # hold beta values over time
n_cells_list = [] # hold # cells over time


total_before = process_time()
i = k = 0    # i is day, k is iteration (i \neq j iff the simulation restarts)
while (not ensmbl.atCapacity()) and (i < total_days+1):
    if k == 1e9:     # don't loop forever
        sys.exit()
    
    print(f'Day {i}')
    
    # Print progress
    if (i > 0) and (i % 50 == 0):
        n_cells = ensmbl.getNumCells()
        print(f'{n_cells} cells')
        print()
    
    beta_list.append(ensmbl.getBetaValues())
    n_cells_list.append(ensmbl.getNumCells())
    
    
    if ensmbl.passDay(): # passDay was successfull
        i += 1
    else:                # had to restart (e.g. all cells died)
        i = 0
        ensmbl.reInit()
        beta_list.clear()
        n_cells_list.clear()

    k += 1
            

after_total = process_time()
print(f'Total time: {after_total - total_before}')

beta_arr = np.stack(beta_list, axis=0)
np.savetxt(os.path.join(prog_params['output_dir'], 'beta_values.txt'), beta_arr, delimiter='\t')
np.savetxt(os.path.join(prog_params['output_dir'], 'n_cells.txt'), n_cells_list, fmt='%d')