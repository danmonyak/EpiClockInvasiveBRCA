import numpy as np
import pandas as pd
import sys
import os
from time import process_time
import simulation_array as sim

output_dir = 'three_sites'

FLIP_RATE = 0.004
GROWTH_RATE = 0.16
DEATH_RATE = 0.14

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
    
    
init_params = {'n_cells':1,
               'flip_rate':FLIP_RATE, # flip rate per cell division per allele
               'growth_rate':GROWTH_RATE, # cell divisions per day
               'death_rate':DEATH_RATE, # cell deaths per day,
               'init_site_state_counts':[1, 1, 0, 1],
              }
gen = np.random.default_rng(0)

nyears = 2

total_days = int(nyears * 365)
ensmbl = sim.Ensemble(init_params, gen)

beta_list = [[], [], []]
n_cells_list = []

# f = open(os.path.join(output_dir, 'progress.txt'), 'w')

total_before = process_time()
i = k = 0
while i < total_days+1:
    if k == 1e9:
        sys.exit()
    
    if (i > 0) and (i % 30 == 0):
        print(i)
        print(ensmbl.getNumCells())
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
np.savetxt(os.path.join(output_dir, 'beta_values.txt'), beta_arr, delimiter='\t')
np.savetxt(os.path.join(output_dir, 'n_cells.txt'), n_cells_list, fmt='%d')