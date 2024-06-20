from time import process_time, sleep
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random

MAX_CELLS = int(1e10)

class Ensemble:
    def __init__(self, init_params, gen):
        self.init_params = init_params
        self.gen = gen
        self.time_obj = dict(zip(['a1', 'a21', 'a22', 'a3', 'b', 'c', 'd'], [0]*8))
        
        self.available_cells = list(range(MAX_CELLS))[::-1]
        self.living_cells = [self.available_cells.pop()]
        
        self.state_arr = np.zeros([MAX_CELLS, sum(self.init_params['init_site_state_counts'])], dtype='byte')
        print('allocated')
        self.reInit()
    def getNumCells(self):
        return len(self.living_cells)
    def reInit(self):
        if ('n_CpGs' in self.init_params) and ('init_site_state_probs' in self.init_params):
#             self.sites = [Site(self.init_params, self.gen) for i in range(self.init_params['n_CpGs'])]
            sys.exit('not set up yet...')
        elif ('init_site_state_counts' in self.init_params) and (self.init_params['n_cells'] == 1):

            
            
            count = 0
            for st in range(4):
                next_count = count + self.init_params['init_site_state_counts'][st]
                self.state_arr[self.living_cells[0], count:next_count] = st
                count = next_count
#             self.state_arr = np.reshape(np.concatenate([[st for i in range(self.init_params['init_site_state_counts'][st])] for st in range(4)], dtype=int, casting='unsafe'), [1, -1])
        else:
            sys.exit('Must provide a init_site_state_counts argument...')
    def getBetaValues(self):
        meth_alleles = np.ceil(self.state_arr[self.living_cells] / 2)
        return (meth_alleles/2).mean(axis=0)
    def passDay(self):
        """
        Cumulative probabilities
        First range - divide
        Second range - die
        Third range - stay
        """
#         event_arr = self.gen.choice(3, size=self.getNumCells(), replace=True, p=[self.init_params['growth_rate'], self.init_params['death_rate'], 1 - self.init_params['growth_rate'] - self.init_params['death_rate']])
#         new_cell_states = self.state_arr[event_arr == 0]
#         self.state_arr = self.state_arr[event_arr != 1, ]
#         if np.all(event_arr == 1):   # tumor was eliminated
#             return False
        
        n_divide, n_die, n_nothing = self.gen.multinomial(self.getNumCells(), pvals=[self.init_params['growth_rate'], self.init_params['death_rate'], 1 - self.init_params['growth_rate'] - self.init_params['death_rate']])
        
#         print(f'n_cells: {self.getNumCells()}')
#         print(f'n_divide: {n_divide}')
#         sleep(1)
        
        if n_divide + n_nothing == 0:   # tumor was eliminated
#             print()
            print('died')
#             print(n_die)
#             print(self.getNumCells())
#             print()

            self.available_cells.extend(self.living_cells)
            self.living_cells = [self.available_cells.pop()]
            return False
        
        self.time_obj['a1'] -= process_time()
        
        self.gen.shuffle(self.living_cells)
        
        self.time_obj['a1'] += process_time()
        self.time_obj['a21'] -= process_time()
        a = self.living_cells[n_divide + n_nothing:]
        
        self.time_obj['a21'] += process_time()
        self.time_obj['a22'] -= process_time()
        
        self.available_cells.extend(a)
        
        self.time_obj['a22'] += process_time()
        self.time_obj['a3'] -= process_time()
        
        self.living_cells = self.living_cells[:n_divide + n_nothing]
        
        self.time_obj['a3'] += process_time()
        self.time_obj['b'] -= process_time()
        
        if n_divide > len(self.available_cells):
            print('Simulation reached capacity')
            self.available_cells.extend(self.living_cells)
            self.living_cells = [self.available_cells.pop()]
            return False
        
        new_cell_states = self.state_arr[self.living_cells[:n_divide]]
        new_cell_idxs = [self.available_cells.pop() for i in range(n_divide)]
        self.living_cells.extend(new_cell_idxs)
        
        
        mu = self.init_params['flip_rate']
        pvals = [(1 - mu)**2, mu*(1 - mu), mu*(1 - mu), mu**2]
        
        self.time_obj['b'] += process_time()
        self.time_obj['c'] -= process_time()
        
        flip_event_arr = self.gen.choice(4, size=new_cell_states.shape, replace=True, p=pvals)
        
        self.time_obj['c'] += process_time()
        self.time_obj['d'] -= process_time()
        
        self.state_arr[new_cell_idxs] = new_cell_states ^ flip_event_arr
#         new_cell_states ^ flip_event_arr
#         self.state_arr = np.vstack([self.state_arr, new_cell_states])
        
        self.time_obj['d'] += process_time()
    
        return True    # tumor is still alive

def plotBetaValues(ax, cell_count_list=None, beta_values=None, binwidth=None, color=None, opacity=None,
                   labelfontsize=None, ticksfontsize=None, sf=1, bins=None):
    if cell_count_list is not None:
        beta_values = [getBetaValues(cell_counts) for cell_counts in cell_count_list]
    elif beta_values is None:
        sys.exit('One must be not None')
    
    sns.histplot(ax=ax, x=beta_values,
#                  binwidth=0.1,
                 bins=bins,
                 stat='proportion', color=color, alpha=opacity)
    ax.set_xlabel('β', fontsize=labelfontsize * sf)
    ax.set_ylabel('Proportion', fontsize=labelfontsize * sf)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlim(-0.05, 1.05)
    ax.tick_params(axis='both', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)
