"""

mRNA
"""
# add plotting scripts' folder as system path
import importlib
import sys
sys.path.append('/data/enguang/CMEODE/analysis/wcm_analyses')

import WCM_analysis

import numpy as np

def plot_average_states_single(w, locusNum, fig_dir, fig_label):
    """
    plot the ensemble averaged count of single mRNA's three states
    """
    prefixes = ['R_', 'Ribosome_mRNA_', 'Degradosome_mRNA_']

    ids = [prefix + locusNum for prefix in prefixes]

    traces = w.get_species_traces(ids)

    title = 'mRNA_{0}_3States'.format(locusNum)

    w.plot_ensemble_averaged_multiples(fig_dir, fig_label, '.png', np.mean(traces, axis=2), ids, 'count', title, True)

    return None

# def plot_hist_state_time(w, locusNum, time_moment, state, fig_dir, fig_label):
#     """
#     plot the distribution of one single mRNA among multiple replicates in certain state at time time_moment
#     time_moment is an integer
#     """
#     state_prefix = {'free':['R_'], 'degradosome':['Degradosome_mRNA_'], 'ribosome':['Ribosome_mRNA_'], 'total':['R_', 'Ribosome_mRNA_']}

#     prefixes = state_prefix[state]

#     value = np.zeros((w.N_reps))

#     for prefix in prefixes: 
#         id = prefix + locusNum
#         value += w.get_specie_trace(id)[time_moment,:]
    
#     xlabel = 'Count'; ylabel = 'Frequency'
#     title = 'Distribution_{0}_mRNA_{1}_time_{2}'.format(state, locusNum, time_moment)

#     w.plot_hist(fig_dir, fig_label, '.png', value, xlabel, ylabel, title, bins=50)


#     return None

def plot_hist_state_time_poisson(w, locusNum, time_moment, state, fig_dir, fig_label):
    """
    plot the distribution of one single mRNA among multiple replicates in certain state at time time_moment
    time_moment is an integer
    """
    state_prefix = {'free':['R_'], 'degradosome':['Degradosome_mRNA_'], 'ribosome':['Ribosome_mRNA_'], 'total':['R_', 'Ribosome_mRNA_']}

    prefixes = state_prefix[state]

    value = np.zeros((w.N_reps))

    for prefix in prefixes: 
        id = prefix + locusNum
        value += w.get_specie_trace(id)[time_moment,:]
    
    xlabel = 'Count'; ylabel = 'Frequency'
    title = 'Distribution_{0}_mRNA_{1}_time_{2}_Poisson'.format(state, locusNum, time_moment)
    w.plot_hist_poisson(fig_dir, fig_label, '.png', value, xlabel, ylabel, title)


    return None

def plot_hist_state_timeperiod_poisson(w, locusNum, time_array, state, fig_dir, fig_label):
    """
    plot the distribution of one single mRNA among multiple replicates in certain state at time time_moment
    time_array is a 1D numpy array
    """
    state_prefix = {'free':['R_'], 'degradosome':['Degradosome_mRNA_'], 'ribosome':['Ribosome_mRNA_'], 'total':['R_', 'Ribosome_mRNA_', 'Degradosome_mRNA_']}

    prefixes = state_prefix[state]

    value = np.zeros((len(time_array), w.N_reps))

    start_time = time_array[0]; end_time = time_array[-1]

    for prefix in prefixes: 
        id = prefix + locusNum
        value += w.get_specie_trace(id)[start_time:end_time+1, :]
    
    value_flattened = value.flatten()

    xlabel = 'Count'; ylabel = 'Frequency'
    title = 'Distribution_{0}_mRNA_{1}_timeperiod_{2}_{3}_Poisson'.format(state, locusNum, start_time, end_time)

    w.plot_hist_poisson(fig_dir, fig_label, '.png', value_flattened, xlabel, ylabel, title)


    return None
