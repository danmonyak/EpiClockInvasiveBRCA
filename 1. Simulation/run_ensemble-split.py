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

def arr2dToString(arr):
    return '\n'.join(['\t'.join(map(lambda x: format(x, '.18e'), x)) for x in arr]) + '\n'

# No index provided
if len(sys.argv) == 1:
    sys.exit('Enter index of parameter dict...')

params_idx = int(sys.argv[1])
params_list = [
    {'output_dir':'90_sites_NB_split', 'n_CpGs_each':30, 'flip_rate':0.005, 'death_rate':0.12, 'nyears':0.5, 'seed':0},
    {'output_dir':'90_sites_2_years_split', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':0},
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
                      split=True, split_limit=int(1e6), n_split=10)

ensmbl_con = sim.EnsembleContainer()
ensmbl_con.addEnsemble(ensmbl)

total_days = int(prog_params['nyears'] * 365)

beta_list = []   # hold beta values over time
n_cells_list = [] # hold # cells over time

beta_values_outfilepath = os.path.join(prog_params['output_dir'], 'beta_values.txt')
n_cells_outfilepath = os.path.join(prog_params['output_dir'], 'n_cells.txt')

if os.path.exists(beta_values_outfilepath):
    os.remove(beta_values_outfilepath)
if os.path.exists(n_cells_outfilepath):
    os.remove(n_cells_outfilepath)

###################################
time_spent = {
    'passDay':0,
    'getBetaValues':0,
    'getBetaValues - line1':0,
    'getBetaValues - line2':0
}
###################################

total_before = process_time()
i = k = 0    # i is day, k is iteration (i \neq j iff the simulation restarts)
found_first = found_second = False
while (not ensmbl_con.atCapacity()) and (i < total_days+1):
    if k == 1e9:     # don't loop forever
        sys.exit()
    
    print(f'Day {i}')
    
    # Print progress
    if (i > 0) and (i % 50 == 0):
        n_cells = ensmbl_con.getNumCells()
        print(f'{n_cells} cells')
        print()
        print(f'Size of state_arr: {ensmbl.state_arr.shape}')

        print('Flushing results to files...')
        
        with open(beta_values_outfilepath, 'a') as f:
            beta_arr = np.stack(beta_list, axis=0)
            beta_arr_str = arr2dToString(beta_arr)
            f.writelines(beta_arr_str)
            beta_list.clear()
        with open(n_cells_outfilepath, 'a') as f:
            n_cells_list_str = '\n'.join(map(str, n_cells_list)) + '\n'
            f.writelines(n_cells_list_str)
            n_cells_list.clear()
    
    time_spent['getBetaValues'] -= time()
    beta_list.append(ensmbl_con.getBetaValues(time_spent))
    time_spent['getBetaValues'] += time()
    
    n_cells_list.append(ensmbl_con.getNumCells())
    
    time_spent['passDay'] -= time()
    response = ensmbl_con.passDay()
    time_spent['passDay'] += time()
    
    if response: # passDay was successful
        i += 1
    else:                # had to restart (e.g. all cells died)
        i = 0
        ensmbl_con.reInit()
        beta_list.clear()
        n_cells_list.clear()

        if os.path.exists(beta_values_outfilepath):
            os.remove(beta_values_outfilepath)
        if os.path.exists(n_cells_outfilepath):
            os.remove(n_cells_outfilepath)

    k += 1



    bound1 = 1e5
    bound2 = 1e6
    n_cells = ensmbl_con.getNumCells()
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

beta_arr_str = arr2dToString(beta_arr)
n_cells_list_str = '\n'.join(map(str, n_cells_list)) + '\n'


# os.remove(beta_values_outfilepath)
# os.remove(n_cells_outfilepath)

# np.savetxt(beta_values_outfilepath, beta_arr, delimiter='\t')
# np.savetxt(n_cells_outfilepath, n_cells_list, fmt='%d')

with open(beta_values_outfilepath, 'a') as f:
    f.writelines(beta_arr_str)

with open(n_cells_outfilepath, 'a') as f:
    f.writelines(n_cells_list_str)
