"""
simulation.py
=======
Author - Daniel Monyak
9-3-24
=======

Source code for simulating fCpG sites in a growing population of cancer cells
Population is modeled as a 2d array of states:
    # cells x # CpGs
    4 possible states (two alleles)
Time is handled discretely
    Pass one day at a time
    Exponential RVs (time to next flip/division/death) are approximated as Bernoulli RVs (probability event occurs in a day)
    Expo(theta) ----> Bernouli(theta)
    
"""

from time import process_time, sleep
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from scipy.stats import linregress

MAX_CELLS = int(1e8)


def betasFromStates(state_arr):
    """
    Returns array of beta values
    
    Parameters
    ----------
    state_arr : 2d array of states of cells x sites
    
    Returns
    -------
    beta_values : ndarray of beta values of fCpGs
    
    """
    meth_alleles = np.ceil(state_arr / 2)
    beta_values = (meth_alleles/2).mean(axis=0)
    return beta_values

def betaCorr(beta1, beta2):
    """
    Returns correlation of beta values
    
    Parameters
    ----------
    beta1, beta2 : arrays of beta values
    
    Returns
    -------
    rvalue : R-value of Pearson Correlation between beta1 and beta2 arrays
    
    """
    rvalue = linregress(beta1, beta2).rvalue
    return rvalue


class Ensemble:
    """
    Contains 2d array of states
    Keeps pointers to
        initial parameters
        random generator object
    Represents an ensemble of fCpG sites in a growing population of cells
    
    Attributes
    ----------
    init_params : dict
        Initial parameters of simulation
    gen : Numpy random Generator
        Used to simulate random events
    max_cells: int
        Maximum number of cells in population
    available_cells: list
        List of cells (represented as ints) that are available for use
    living_cells: list
        List of living cells
    state_arr: 2d ndarray
        # cells x # CpGs
        4 possible states
            0: both alleles unmethylated
            1: first allele methylated, second allele unmethylated
            2: first allele unmethylated, second allele methylated
            3: both alleles methylated
            
            These state definitions derive their binary representations
                0 = 00b
                1 = 01b
                2 = 10b
                3 = 11b
                
                0 = unmethylated, 1 = methylated
                
    at_capacity: bool
        All possible cells have been exhausted
        
    Methods
    -------
    
    """
    
    def __init__(self, init_params, gen, max_cells=MAX_CELLS):
        """
        Parameters
        ----------
        init_params : dict
            Initial parameters of simulation
        gen : Numpy random Generator
            Used to simulate random events
        max_cells: int
            Maximum number of cells in population
            
        """
        
        self.init_params = init_params
        self.gen = gen
        self.max_cells = max_cells
        
        self.available_cells = list(range(self.max_cells))[::-1]
        self.living_cells = [self.available_cells.pop()]        # Initialize with one cell
        
        # Determine # CpG sites
        if 'n_CpGs' in self.init_params:
            n_CpGs = self.init_params['n_CpGs']
        else:
            n_CpGs = sum(self.init_params['init_site_state_counts'])
        
        
        self.state_arr = np.zeros([self.max_cells, n_CpGs], dtype='byte')
        print('Allocated memory!')
        
        self.reInit()
        self.at_capacity = False
    
    def copy(self, new_gen):
        """
        Create copy of fCpG ensemble
        Requires passing of a new generator object, so that copied ensembles do not generate the same random events
        
        Parameters
        ----------
        new_gen : Numpy random Generator
            New generator object
            
        Returns
        -------
        new_ens : Ensemble
            Copy of this Ensemble object
            
        """
        
        new_ens = Ensemble(self.init_params, new_gen, self.max_cells)
        new_ens.available_cells = self.available_cells.copy()
        new_ens.living_cells = self.living_cells.copy()
        new_ens.state_arr = self.state_arr.copy()
        return new_ens
    
    def getRandCell(self, max_cells=None):
        """
        Create Ensemble representing a metastasis of the current population of cells
        Create a new Ensemble and copy the state of a randomly picked living cell
        
        Parameters
        ----------
        max_cells : int
            maximum # cells in new Ensemble
            
        Returns
        -------
        new_ens : Ensemble
            New metastasis Ensemble
            
        """
        if max_cells is None:
            max_cells = self.max_cells
        new_ens = Ensemble(self.init_params, np.random.default_rng(), max_cells)
        cell_i = np.random.choice(self.living_cells, replace=False)
        new_ens.state_arr[new_ens.living_cells[0]] = self.state_arr[cell_i]
        return new_ens
    
    def getNumCells(self):
        """
        Return the number of living cells
        """
        return len(self.living_cells)
    
    def atCapacity(self):
        """
        Return True iff Ensemble has reached max number of living cells
        """
        return self.at_capacity
    
    def reInit(self):
        """
        Initialize or reinitialize (if all cells died) Ensemble 
        Only initialize one cells (that's how many are alive)
        
        If n_CpGs and init_site_state_probs are both in self.init_params
            initialize first cell with that number of CpGs and randomly distribute it among
            the 4 states according to the probabilities in init_site_state_probs
        If not, but init_site_state_counts is in self.init_params
            init_site_state_counts holds desired number of sites in each state in first cell
        """
        
        if ('n_CpGs' in self.init_params) and ('init_site_state_probs' in self.init_params):
            self.state_arr[self.living_cells[0]] = self.gen.choice(4, p=self.init_params['init_site_state_probs'], size=[1, self.init_params['n_CpGs']])
        elif 'init_site_state_counts' in self.init_params:
            count = 0
            for st in range(4):
                next_count = count + self.init_params['init_site_state_counts'][st]
                self.state_arr[self.living_cells[0], count:next_count] = st
                count = next_count
#             self.state_arr = np.reshape(np.concatenate([[st for i in range(self.init_params['init_site_state_counts'][st])] for st in range(4)], dtype=int, casting='unsafe'), [1, -1])
        else:
            sys.exit('Must provide a init_site_state_counts argument...')
    
    def getBetaValues(self):
        """
        Return beta values of all fCpGs, using all living cells
            
        Returns
        -------
        beta_values : ndarray
            Beta values of all fCpG sites
            
        """
        beta_values = betasFromStates(self.state_arr[self.living_cells])
        return beta_values
    
    def getRandCellsBetaValues(self, n=1):
        """
        Return beta values of n randomly picked living cells
            
        Returns
        -------
        beta_value_list : List of ndarrays
            Each ndarray is the array of beta values of a single cell
            
        """
        
        n = min(n, self.getNumCells())
        beta_value_list = []
        for cell in np.random.choice(self.living_cells, size=n, replace=False):
            beta_value_list.append(betasFromStates(np.reshape(self.state_arr[cell], [1, -1])))
        
        return beta_value_list
            
    def passDay(self, die=False):
        """
        Pass one time unit (one day)
        
        Pick number of each event from multinomial model
            Cumulative probabilities for events
            First range - divide
            Second range - die
            Third range - nothing
        If tumor is eliminated, reset the Ensemble and return False
        Dying cells - return cells to self.available_cells
        Dividing cells
            - pick cells from self.available_cells and copy states
            - determine new state of each new cell using flip rate to determine probabilities
        
        Parameters
        ----------
        die : bool
            if True, it forces the Ensemble has to die
        
        Returns
        -------
           bool : False if the tumor was eliminated or it ran out of available cells, else True
        
        """
        
        n_divide, n_die, n_nothing = self.gen.multinomial(self.getNumCells(), pvals=[self.init_params['growth_rate'], self.init_params['death_rate'], 1 - self.init_params['growth_rate'] - self.init_params['death_rate']])
        
        if die or (n_divide + n_nothing == 0):   # Tumor was eliminated
            print('Tumor died')
            self.available_cells.extend(self.living_cells)
            self.living_cells = [self.available_cells.pop()]
            return False
        
        # Shuffle order of living cells so that the same cells don't live forever
        self.gen.shuffle(self.living_cells)
        
        # Return dead cells to self.available_cells
        self.available_cells.extend(self.living_cells[n_divide + n_nothing:])
        self.living_cells = self.living_cells[:n_divide + n_nothing]
        
        # If the Ensemble will run out of cells
        # Reset cell lists and return False
        if n_divide > len(self.available_cells):
            print('Simulation reached capacity')
            self.at_capacity = True
            self.available_cells.extend(self.living_cells)
            self.living_cells = [self.available_cells.pop()]
            return False
        
        # Create new cells from divisions
        new_cell_states = self.state_arr[self.living_cells[:n_divide]]
        new_cell_idxs = [self.available_cells.pop() for i in range(n_divide)]
        self.living_cells.extend(new_cell_idxs)
        
        # Probability of each event
        #       Neither allele flips, first allele flips, second allele flips, both alleles flip
        mu = self.init_params['flip_rate']
        pvals = [(1 - mu)**2, mu*(1 - mu), mu*(1 - mu), mu**2]
        
        # Holds each state change for each new cell
        flip_event_arr = self.gen.choice(4, size=new_cell_states.shape, replace=True, p=pvals)
        
        # Bitwse XOR operation correctly performs state transitions
        #      0 (unmethylated) ^ 0 (no flip) = 0 (unmethylated)
        #      0 (unmethylated) ^ 1 (flip) = 1 (methylated)
        #      1 (methylated) ^ 0 (no flip) = 1 (methylated)
        #      1 (methylated) ^ 1 (flip) = 0 (unmethylated)
        self.state_arr[new_cell_idxs] = new_cell_states ^ flip_event_arr
#         new_cell_states ^ flip_event_arr
#         self.state_arr = np.vstack([self.state_arr, new_cell_states])
    
        return True    # tumor is still alive


def plotBetaValues(ax, beta_values=None, binwidth=None, color=None, opacity=None,
                   labelfontsize=None, ticksfontsize=None, sf=1, bins=None):
    """
    Plot histogram of beta values

    Parameters
    ----------
    ax : Matplotlib Axes object
        Axes to be used for plot
    beta_values : ndarray
        Beta values to be plotted
    binwidth : number or pair of numbers
    color : str or tuple of floats
    opacity : float
    labelfontsize : float
    ticksfontsize : float
    sf : float
        scale factor of figure
    bins : str, number, vector, or a pair of such values
        See Seaborn histplot documentation

    Returns
    -------
       bool : False if the tumor was eliminated or it ran out of available cells, else True

    """
    
    sns.histplot(ax=ax, x=beta_values, bins=bins, stat='proportion', color=color, alpha=opacity)
    ax.set_xlabel('β', fontsize=labelfontsize * sf)
    ax.set_ylabel('Proportion', fontsize=labelfontsize * sf)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlim(-0.05, 1.05)
    ax.tick_params(axis='both', labelsize=ticksfontsize * sf, width=sf, length=8 * sf)
