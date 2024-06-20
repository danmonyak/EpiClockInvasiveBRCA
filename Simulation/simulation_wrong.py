import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random

def threeSidedFlip(pvals):
    r = random.random()
    if r > (pvals[0] + pvals[1]):
        return 2
    if r > pvals[0]:
        return 1
    return 0

def pickEvents(n, event_counts, gen):
    ret = np.zeros(3, dtype=int)
    for i in range(n):
        try:
#             choice = gen.choice(3, p=event_counts / event_counts.sum())
            choice = threeSidedFlip(event_counts / event_counts.sum())
        except:
            print(event_counts)
            sys.exit()
            
        event_counts[choice] -= 1
        ret[choice] += 1
    
    return ret

class Site:
    def __init__(self, init_params, gen, state=None):
        self.init_params = init_params
        self.gen = gen
        if state is None:
            self.initCells()
        else:
            self.initOneCell(state)
    def initCells(self):
        self.cell_counts = np.zeros(4, dtype=int)
        for state in self.gen.choice(4, p=self.init_params['init_site_state_probs'], size=self.init_params['n_cells']):
            self.cell_counts[state] += 1
    def initOneCell(self, state):
        self.cell_counts = np.zeros(4, dtype=int)
        self.cell_counts[state] += 1
    def getBetaValue(self):
        n_cells = np.sum(self.cell_counts)
        meth_alleles = self.cell_counts[1] + self.cell_counts[2] + 2 * self.cell_counts[3]
        return meth_alleles / (2 * n_cells)
    def passDay(self, event_arr):
#     def passDay(self, event_counts):
        cell_counts_delta = [0]*4

        prev_count = 0
        for state in range(4):
            next_count = prev_count + self.cell_counts[state]
            event_arr_state = event_arr[prev_count : next_count]
            prev_count = next_count
            n_divide = np.sum(event_arr_state == 0)
            n_die = np.sum(event_arr_state == 1)
#             n_divide, n_die, _ = pickEvents(self.cell_counts[state], event_counts, self.gen)
    
            cell_counts_delta[state] -= n_die

            mu = self.init_params['flip_rate']
            pvals=[mu**2, mu*(1 - mu), mu*(1 - mu), (1 - mu)**2]
            assert abs(np.sum(pvals) - 1) < 1e-07
            flip_both, flip_allele_one, flip_allele_two, flip_none = self.gen.multinomial(n_divide, pvals=pvals)

            cell_counts_delta[state ^ 3] += flip_both
            cell_counts_delta[state ^ 1] += flip_allele_one
            cell_counts_delta[state ^ 2] += flip_allele_two
            cell_counts_delta[state] += flip_none
            
#             if (flip_both + flip_allele_one + flip_allele_two) > 0:
#                 print('flip')

        self.cell_counts += cell_counts_delta

class Ensemble:
    def __init__(self, init_params, gen):
        self.init_params = init_params
        self.gen = gen
        self.n_cells = init_params['n_cells']
        self.reInit()
    def reInit(self):
        if ('n_CpGs' in self.init_params) and ('init_site_state_probs' in self.init_params):
            self.sites = [Site(self.init_params, self.gen) for i in range(self.init_params['n_CpGs'])]
        elif ('init_site_state_counts' in self.init_params) and (self.init_params['n_cells'] == 1):
            self.sites = []
            for state in range(4):
                self.sites.extend([Site(self.init_params, self.gen, state) for i in range(self.init_params['init_site_state_counts'][state])])
        else:
            sys.exit('Must provide a init_site_state_counts argument...')
    def getBetaValues(self):
        return np.array([s.getBetaValue() for s in self.sites])
    def passDay(self):
        """
        Cumulative probabilities
        First range - divide
        Second range - die
        Third range - stay
        """
        event_arr = self.gen.choice(3, size=self.n_cells, replace=True, p=[self.init_params['growth_rate'], self.init_params['death_rate'], 1 - self.init_params['growth_rate'] - self.init_params['death_rate']])
        n_divide = np.sum(event_arr == 0)
        n_die = np.sum(event_arr == 1)
#         event_counts = self.gen.multinomial(self.n_cells, pvals=[self.init_params['growth_rate'], self.init_params['death_rate'], 1 - self.init_params['growth_rate'] - self.init_params['death_rate']])
#         n_divide, n_die, _ = event_counts
        
        if n_die == self.n_cells:   # tumor was eliminated
            return False  
        
        self.n_cells += n_divide - n_die
        
        for s in self.sites:
            s.passDay(event_arr)
#             s.passDay(event_counts.copy())
    
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
