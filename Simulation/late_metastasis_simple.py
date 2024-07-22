import numpy as np
import pandas as pd
import sys
import os
from time import process_time
import simulation_array as sim

FLIP_RATE = 0.004
GROWTH_RATE = 0.16
DEATH_RATE = 0.14

output_dir = 'late_met'

if not os.path.exists(output_dir):
    os.mkdir(output_dir)
    
    
init_params = {'flip_rate':FLIP_RATE, # flip rate per cell division per allele
               'growth_rate':GROWTH_RATE, # cell divisions per day
               'death_rate':DEATH_RATE, # cell deaths per day,
               'n_CpGs':90,
               'init_site_state_probs':[1/3, 1/3, 0, 1/3],
              }

nyears = 2
total_days = int(nyears * 365)

ensmbl = sim.Ensemble(init_params, np.random.default_rng(0))

beta_list = []
n_cells_list = []
years_list = []

# beta_one_cell_list = []

rvalue_list = []
num_corr = 30

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
    
    bulk_beta = ensmbl.getBetaValues()
    n_cells = ensmbl.getNumCells()
    if (i % 10 == 0) and (n_cells >= num_corr):
        n_cells_list.append(n_cells)
        years_list.append(i / 365)
        rvalue_list.append([sim.betaCorr(bulk_beta, single_cell_beta) for single_cell_beta in ensmbl.getRandCellsBetaValues(num_corr)])

    before = process_time()

    if ensmbl.passDay():
        i += 1
    else:
        i = 0
        ensmbl.reInit()
        beta_list.clear()
        n_cells_list.clear()
        rvalue_list.clear()
        years_list.clear()

    after = process_time()

    day_dur = after - before
    total_dur = after - total_before

    k += 1

after_total = process_time()
print(f'Total time: {after_total - total_before}')

# beta_arr = np.stack(beta_list, axis=0)
# beta_one_cell_arr = np.stack(beta_one_cell_list, axis=0)

# np.savetxt(os.path.join(output_dir, 'beta_values.txt'), beta_arr, delimiter='\t')
# np.savetxt(os.path.join(output_dir, 'beta_one_cell.txt'), beta_one_cell_arr, fmt='%d')

# r_values_df = pd.DataFrame(index=n_cells_list, data=np.array(rvalue_list), columns=[f'cell_{i}' for i in range(num_corr)])
# r_values_df.index.name = 'n_cells'

# r_values_df = pd.concat([pd.Series(years_list, name='year'), pd.DataFrame(np.array(rvalue_list), columns=[f'cell_{i}' for i in range(num_corr)])], axis=1)

r_values_df = pd.DataFrame(np.array(rvalue_list).T)

r_values_df.to_csv(os.path.join(output_dir, 'r_values.txt'), sep='\t')
np.savetxt(os.path.join(output_dir, 'years_list.txt'), years_list, fmt='%f')
