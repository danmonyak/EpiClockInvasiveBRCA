import numpy as np
import pandas as pd
import sys
import os
from time import process_time
import simulation_array as sim

FLIP_RATE = 0.004
GROWTH_RATE = 0.16
DEATH_RATE = 0.14

output_dir = 'synchro_met'

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# FLIP_RATE /= 2
# DEATH_RATE += 0.006
    
init_params = {'flip_rate':FLIP_RATE, # flip rate per cell division per allele
               'growth_rate':GROWTH_RATE, # cell divisions per day
               'death_rate':DEATH_RATE, # cell deaths per day,
               'n_CpGs':90,
               'init_site_state_probs':[1/3, 1/3, 0, 1/3],
              }

nyears = 2
total_days = int(nyears * 365)

ensmbl_1 = sim.Ensemble(init_params, np.random.default_rng(0), max_cells=int(1e7))
ensmbl_1_met = ensmbl_1.copy(np.random.default_rng(1))
ensmbl_2 = sim.Ensemble(init_params, np.random.default_rng(2), max_cells=int(1e7))

# ensmbl_1 = sim.Ensemble(init_params, np.random.default_rng(0), max_cells=int(1e7))
# ensmbl_1_met = ensmbl_1.copy(np.random.default_rng(1))

init_state = ensmbl_1.state_arr[ensmbl_1.living_cells[0]].copy()

total_before = process_time()
for label, ensmbl in [('ensmbl_1', ensmbl_1), ('ensmbl_1_met', ensmbl_1_met),
#                       ('ensmbl_2', ensmbl_2)
                     ]:
    beta_list = []
    n_cells_list = []

    i = k = 0
    while (not ensmbl.atCapacity()) and (i < total_days+1):
        n_cells = ensmbl.getNumCells()
        
        if k == 1e9:
            sys.exit()

        print(f'{label}, day {i}')
        if (i > 0) and (i % 50 == 0):
            print(f'{n_cells} cells')
            print()

        beta_list.append(ensmbl.getBetaValues())
        n_cells_list.append(n_cells)

        before = process_time()

        if ensmbl.passDay():
            i += 1
        else:
            i = 0
            ensmbl.reInit()
            ensmbl.state_arr[ensmbl.living_cells[0]] = init_state
            beta_list.clear()
            n_cells_list.clear()

        after = process_time()

        day_dur = after - before
        total_dur = after - total_before

        k += 1

    del ensmbl
    after_total = process_time()
    print(f'Total time: {after_total - total_before}')

    beta_arr = np.vstack(beta_list)
    np.savetxt(os.path.join(output_dir, f'beta_values_{label}.txt'), beta_arr, delimiter='\t')
    np.savetxt(os.path.join(output_dir, f'n_cells_{label}.txt'), n_cells_list, fmt='%d')