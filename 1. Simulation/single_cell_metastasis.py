"""
single_cell_metastasis.py
=======
Author - Daniel Monyak
9-5-24
=======

Source code for running a simulation of a simplified version of a late metastasis using an fCpG Ensemble
A metastasis in this case is just a single cell chosen from the primary tumor
At each time point, randomly pick a cell and measure Pearson R correlation of its beta values with beta values of the whole tumor

Do this 30 times and record the R value at each time point in each simulatino
"""

import numpy as np
import pandas as pd
import sys
import os
from time import process_time, sleep
import simulation as sim

# Parameters
FLIP_RATE = 0.002
PROLIF_RATE = 0.17
DEATH_RATE = 0.15

output_dir = 'sc_met'

os.makedirs(output_dir, exist_ok=True)
    
init_params = {'flip_rate':FLIP_RATE, # flip rate per cell division per allele
               'prolif_rate':PROLIF_RATE, # cell divisions per day
               'death_rate':DEATH_RATE, # cell deaths per day,
               'n_CpGs':90,
               'init_site_state_probs':[1/3, 1/3, 0, 1/3],
              }

nyears = 2
total_days = int(nyears * 365)

n_sims = 30

# Divide total time into 8 checkpoints
time_checkpoints_days, time_offset_days = np.linspace(0, total_days, 9, retstep=True)
time_checkpoints_days = [int(x) for x in time_checkpoints_days]
time_offset_days = int(time_offset_days)

# Holds calculated rvalue at each time
rvalue_list_list = []

# Each simulation
for sim_i in range(n_sims):
    print('###################')
    print(f'#####  Sim {sim_i}  #####')
    print('###################')
    
    # Create new ensemble
    ensmbl = sim.Ensemble(init_params, np.random.default_rng(), max_cells=int(1e7))
    rvalue_list = []

    total_before = process_time()
    i = k = 0   # i is day, k is iteration (i \neq j iff the simulation restarts)

    while i < total_days+1:
        if k == 1e9:
            sys.exit()

        print(f'Sim {sim_i}, Day {i}')
        
        # Print progress
        if (i > 0) and (i % 50 == 0):
            n_cells = ensmbl.getNumCells()
            print(f'{n_cells} cells')
            print()
        
        # Time checkpoint
        # Calculate and record R value between beta values of single cell and entire tumor
        if i % time_offset_days == 0:
            bulk_beta = ensmbl.getBetaValues()
            single_cell_beta = ensmbl.getRandCellsBetaValues()
            rvalue_list.append(sim.betaCorr(bulk_beta, single_cell_beta))
        
        # Tests ensure that simulation will not run too long
        # Necessary assuming compute power
        test = process_time() - total_before > 240
        test = test or ( (process_time() - total_before > 60/10) and (i < 600) )
        test = test or ((i == 100) and (ensmbl.getNumCells() > 20))
        
        if ensmbl.passDay(test):                # passDay was successfull
            i += 1
        else:                                   # had to restart (e.g. all cells died)
            total_before = process_time()
            i = 0
            ensmbl.reInit()
            rvalue_list.clear()
            continue
        
        k += 1

    rvalue_list_list.append(rvalue_list)
    after_total = process_time()
    print(f'Total time: {after_total - total_before}')

r_values_df = pd.DataFrame(index=time_checkpoints_days, data=np.array(rvalue_list_list).T)
r_values_df.to_csv(os.path.join(output_dir, 'r_values.txt'), sep='\t')
