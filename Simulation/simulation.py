import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

N_FLIPS = 0

class Site:
    def __init__(self, init_params, gen):
        self.init_params = init_params
        self.gen = gen
        self.reInit()
    def reInit(self):
        self.cell_counts = initCells(self.init_params, self.gen)

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
    
def getBetaValue(cell_counts):
    n_cells = np.sum(cell_counts)
    meth_alleles = cell_counts[1] + cell_counts[2] + 2 * cell_counts[3]
    return meth_alleles / (2 * n_cells)

def initCells(init_params, gen):
    cell_counts = [0]*4
    for state in gen.choice(4, p=init_params['init_site_state_probs'], size=init_params['n_cells']):
        cell_counts[state] += 1
    return np.array(cell_counts)

def passDay(cell_counts, init_params, gen):
    global N_FLIPS
    
    """
    Cumulative probabilities
    First range - divide
    Second range - die
    Third range - stay
    """
    
    cell_counts_delta = [0]*4
    
    for state in range(4):

        n_divide, n_die, _ = gen.multinomial(cell_counts[state], pvals=[init_params['growth_rate'], init_params['death_rate'], 1 - init_params['growth_rate'] - init_params['death_rate']])
        
        # Prevent cell population from being eliminated at the beginning by chance
        if (cell_counts[state] > 0) and (n_die == cell_counts[state]):
            return False  # tumor was eliminated
        
        cell_counts_delta[state] -= n_die
        
        mu = init_params['flip_rate']
        pvals=[mu**2, mu*(1 - mu), mu*(1 - mu), (1 - mu)**2]
        assert abs(np.sum(pvals) - 1) < 1e-07
        flip_both, flip_allele_one, flip_allele_two, flip_none = gen.multinomial(n_divide, pvals=pvals)
        
        n_flips_raw = (2 * flip_both + flip_allele_one + flip_allele_two)
        if n_flips_raw > 0:
            N_FLIPS += n_flips_raw / cell_counts[state]
        
        cell_counts_delta[state ^ 3] += flip_both
        cell_counts_delta[state ^ 1] += flip_allele_one
        cell_counts_delta[state ^ 2] += flip_allele_two
        cell_counts_delta[state] += flip_none

    cell_counts += cell_counts_delta
    return True    # tumor is still alive