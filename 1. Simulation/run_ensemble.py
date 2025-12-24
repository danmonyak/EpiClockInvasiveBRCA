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
from time import process_time, time
import simulation as sim

# No index provided
if len(sys.argv) == 1:
    sys.exit('Enter index of parameter dict...')

params_idx = int(sys.argv[1])
params_list = [
    {'output_dir':'3_sites', 'n_CpGs_each':1, 'death_rate':0.15, 'nyears':2, 'seed':0},
    {'output_dir':'90_sites', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':0},
    {'output_dir':'larger_death_rate', 'n_CpGs_each':1, 'death_rate':0.16, 'nyears':4, 'seed':0},
    # For neuroblastoma
    {'output_dir':'90_sites_3_years', 'n_CpGs_each':30, 'death_rate':0.16, 'nyears':3, 'seed':0},
    {'output_dir':'90_sites_3_years_2', 'n_CpGs_each':30, 'flip_rate':0.005, 'death_rate':0.14, 'nyears':1.25, 'seed':0},
    ### New simulations
    # delete
    # {'output_dir':'90_sites_seed1', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':1},
    # {'output_dir':'90_sites_seed2', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':2},
    # {'output_dir':'90_sites_seed3', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':3},
    {'output_dir':'new_test_dec23', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':0},
    {'output_dir':'90_sites_NB', 'n_CpGs_each':30, 'flip_rate':0.005, 'death_rate':0.14, 'nyears':1.25, 'seed':0},
]
prog_params = params_list[params_idx]

print(f'Running simulation with the following parameters: {prog_params}')

# Constant parameters
FLIP_RATE = 0.002
PROLIF_RATE = 0.17

if 'flip_rate' not in prog_params:
    prog_params['flip_rate'] = FLIP_RATE
if 'prolif_rate' not in prog_params:
    prog_params['prolif_rate'] = PROLIF_RATE

os.makedirs(prog_params['output_dir'], exist_ok=True)

init_params = {'flip_rate':prog_params['flip_rate'], # flip rate per cell division per allele
               'prolif_rate':prog_params['prolif_rate'], # cell divisions per day
               'death_rate':prog_params['death_rate'], # cell deaths per day,
               'init_site_state_counts':[prog_params['n_CpGs_each'], prog_params['n_CpGs_each'], 0, prog_params['n_CpGs_each']],
              }
gen = np.random.default_rng(prog_params['seed'])
ensmbl = sim.Ensemble(init_params, gen,
                     max_cells=int(2e6))

total_days = int(prog_params['nyears'] * 365)

beta_list = []   # hold beta values over time
n_cells_list = [] # hold # cells over time

###################################
time_spent = {
    'passDay':0,
    'getBetaValues':0,
}
###################################

total_before = process_time()
i = k = 0    # i is day, k is iteration (i \neq j iff the simulation restarts)
found_first = found_second = False
while (not ensmbl.atCapacity()) and (i < total_days+1):
    if k == 1e9:     # don't loop forever
        sys.exit()
    
    print(f'Day {i}')
    
    # Print progress
    if (i > 0) and (i % 50 == 0):
        n_cells = ensmbl.getNumCells()
        print(f'{n_cells} cells')
        print()
        print(f'Size of state_arr: {ensmbl.state_arr.shape}')

    
    time_spent['getBetaValues'] -= time()
    beta_list.append(ensmbl.getBetaValues())
    time_spent['getBetaValues'] += time()
    
    n_cells_list.append(ensmbl.getNumCells())

    time_spent['passDay'] -= time()
    response = ensmbl.passDay()
    time_spent['passDay'] += time()
    
    if response and (response['result'] == 'success'): # passDay was successful
        i += 1
    else:                # had to restart (e.g. all cells died)
        i = 0
        ensmbl.reInit()
        beta_list.clear()
        n_cells_list.clear()

    k += 1

    bound1 = 1e5
    bound2 = 1e6
    n_cells = ensmbl.getNumCells()
    if (not found_first) and (n_cells >= bound1):
        found_first = True
        first_time = time()
    if (not found_second) and (n_cells >= bound2):
        found_second = True
        second_time = time()
        print(f'Took {second_time - first_time} seconds to get from {1e5:,.0f} to {1e6:,.0f} cells')
            

after_total = process_time()
print(f'Total time: {after_total - total_before}')
print(f'Specific time: {time_spent}')

beta_arr = np.stack(beta_list, axis=0)
np.savetxt(os.path.join(prog_params['output_dir'], 'beta_values.txt'), beta_arr, delimiter='\t')
np.savetxt(os.path.join(prog_params['output_dir'], 'n_cells.txt'), n_cells_list, fmt='%d')