import numpy as np
import pandas as pd
import sys
import os
from time import process_time
import simulation_array as sim

if len(sys.argv) == 1:
    sys.exit('Enter index of parameter dict...')

params_idx = int(sys.argv[1])

params_list = [
    {'output_dir':'3_sites', 'n_CpGs_each':1, 'death_rate':0.15, 'nyears':2, 'seed':0},
    {'output_dir':'90_sites', 'n_CpGs_each':30, 'death_rate':0.15, 'nyears':2, 'seed':0},
    {'output_dir':'larger_death_rate', 'n_CpGs_each':1, 'death_rate':0.16, 'nyears':4, 'seed':0},
]

prog_params = params_list[params_idx]

FLIP_RATE = 0.002
GROWTH_RATE = 0.17
# DEATH_RATE = 0.15

os.makedirs(prog_params['output_dir'], exist_ok=True)    
    
init_params = {'flip_rate':FLIP_RATE, # flip rate per cell division per allele
               'growth_rate':GROWTH_RATE, # cell divisions per day
               'death_rate':prog_params['death_rate'], # cell deaths per day,
               'init_site_state_counts':[prog_params['n_CpGs_each'], prog_params['n_CpGs_each'], 0, prog_params['n_CpGs_each']],
              }
gen = np.random.default_rng(prog_params['seed'])

total_days = int(prog_params['nyears'] * 365)
ensmbl = sim.Ensemble(init_params, gen)

beta_list = []
n_cells_list = []

# f = open(os.path.join(prog_params['output_dir'], 'progress.txt'), 'w')

total_before = process_time()
i = k = 0
while (not ensmbl.atCapacity()) and (i < total_days+1):
    if k == 1e9:
        sys.exit()
    
    print(f'Day {i}')
    if (i > 0) and (i % 50 == 0):
        n_cells = ensmbl.getNumCells()
        print(f'{n_cells} cells')
        print()
    
    beta_list.append(ensmbl.getBetaValues())
    n_cells_list.append(ensmbl.getNumCells())

    before = process_time()
    
    if ensmbl.passDay():
        i += 1
    else:
        i = 0
        ensmbl.reInit()
        beta_list.clear()
        n_cells_list.clear()

    after = process_time()
    
    day_dur = after - before
    total_dur = after - total_before
    
#     f.write(f'{i} {ensmbl.getNumCells()} {round(day_dur, 3)} {round(total_dur, 3)}\n')
#     f.flush()
    
    k += 1
            
# f.close()

print(ensmbl.time_obj)

after_total = process_time()
print(f'Total time: {after_total - total_before}')

beta_arr = np.stack(beta_list, axis=0)
np.savetxt(os.path.join(prog_params['output_dir'], 'beta_values.txt'), beta_arr, delimiter='\t')
np.savetxt(os.path.join(prog_params['output_dir'], 'n_cells.txt'), n_cells_list, fmt='%d')