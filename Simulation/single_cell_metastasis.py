import numpy as np
import pandas as pd
import sys
import os
from time import process_time, sleep
import simulation_array as sim

FLIP_RATE = 0.002
GROWTH_RATE = 0.17
DEATH_RATE = 0.15

output_dir = 'sc_met'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
    
init_params = {'flip_rate':FLIP_RATE, # flip rate per cell division per allele
               'growth_rate':GROWTH_RATE, # cell divisions per day
               'death_rate':DEATH_RATE, # cell deaths per day,
               'n_CpGs':90,
               'init_site_state_probs':[1/3, 1/3, 0, 1/3],
              }

nyears = 2
total_days = int(nyears * 365)

n_sims = 30

time_checkpoints_days, time_offset_days = np.linspace(0, total_days, 9, retstep=True)
time_checkpoints_days = [int(x) for x in time_checkpoints_days]
time_offset_days = int(time_offset_days)

# time_checkpoints_days = [int(x * 365) for x in time_checkpoints]
# time_offset_days = int(time_offset * 365)


rvalue_list_list = []

for sim_i in range(n_sims):
    print('###################')
    print(f'#####  Sim {sim_i}  #####')
    print('###################')
#     sleep(2)
    
    ensmbl = sim.Ensemble(init_params, np.random.default_rng(), max_cells=int(1e7))
    rvalue_list = []

    total_before = process_time()
    i = k = 0
#     while (not ensmbl.atCapacity()) and (i < total_days+1):
    while i < total_days+1:
        if k == 1e9:
            sys.exit()

        print(f'Sim {sim_i}, Day {i}')
        if (i > 0) and (i % 50 == 0):
            n_cells = ensmbl.getNumCells()
            print(f'{n_cells} cells')
            print()
            
        if i % time_offset_days == 0:
            bulk_beta = ensmbl.getBetaValues()
            single_cell_beta = ensmbl.getRandCellsBetaValues()
            rvalue_list.append(sim.betaCorr(bulk_beta, single_cell_beta))
        
        test = process_time() - total_before > 240
        test = test or ( (process_time() - total_before > 60/10) and (i < 600) )
        test = test or ((i == 100) and (ensmbl.getNumCells() > 20))
        
        if ensmbl.passDay(test):
            i += 1
        else:
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
