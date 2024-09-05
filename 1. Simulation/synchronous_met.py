"""
synchronous_met.py
=======
Author - Daniel Monyak
9-5-24
=======

Source code for running a simulation of a simplified version of a synchronous metastasis using an fCpG Ensemble
A metastasis in this case is just a second tumor growing synchronously with the primary tumor that
     has the same founding cell
At each time point, randomly pick a cell and measure Pearson R correlation of the beta values of the two tumors

Do this 30 times and record the R value at each time point in each simulatino
"""

import numpy as np
import sys
import os
from time import process_time
import simulation as sim

FLIP_RATE = 0.004
GROWTH_RATE = 0.16
DEATH_RATE = 0.14

output_dir = 'synchro_met'

os.makedirs(output_dir, exist_ok=True)

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
# ensmbl_2 = sim.Ensemble(init_params, np.random.default_rng(2), max_cells=int(1e7))

init_state = ensmbl_1.state_arr[ensmbl_1.living_cells[0]].copy()

total_before = process_time()

# Each tumor
for label, ensmbl in [('ensmbl_1', ensmbl_1), ('ensmbl_1_met', ensmbl_1_met),
#                       ('ensmbl_2', ensmbl_2)
                     ]:
    
    beta_list = []     # hold beta values over time
    n_cells_list = []  # hold # cells over time

    i = k = 0    # i is day, k is iteration (i \neq j iff the simulation restarts)
    while (not ensmbl.atCapacity()) and (i < total_days+1):
        n_cells = ensmbl.getNumCells()
        
        if k == 1e9:
            sys.exit()

        # Print progress
        print(f'{label}, day {i}')
        if (i > 0) and (i % 50 == 0):
            print(f'{n_cells} cells')
            print()

        beta_list.append(ensmbl.getBetaValues())
        n_cells_list.append(n_cells)

        if ensmbl.passDay():                # passDay was successfull
            i += 1
        else:                               # had to restart (e.g. all cells died)
            i = 0
            ensmbl.reInit()
            ensmbl.state_arr[ensmbl.living_cells[0]] = init_state    # need to reset founding cell state
            beta_list.clear()
            n_cells_list.clear()


        k += 1

    del ensmbl
    after_total = process_time()
    print(f'Total time: {after_total - total_before}')

    beta_arr = np.vstack(beta_list)
    np.savetxt(os.path.join(output_dir, f'beta_values_{label}.txt'), beta_arr, delimiter='\t')
    np.savetxt(os.path.join(output_dir, f'n_cells_{label}.txt'), n_cells_list, fmt='%d')