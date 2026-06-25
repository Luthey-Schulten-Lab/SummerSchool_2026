"""
Author: Enguang Fu

Date: March 2024

export the time traces of counts, Surface area, and fluxes into CSV files
"""

import os
import glob

import pandas as pd
import numpy as np

import  integrate
import rxns_ODE as ODE

from math import floor
from math import log10


# =====================================================================================
# Append-only streaming CSV output (single CME run, no restart).
#
# The legacy writeCountstoCSV/SA/Flux re-read the whole growing CSV every restart
# (O(N^2)) and index counts by ABSOLUTE frame number (needs full history). For the
# no-restart design we instead, every output interval, write only the NEW frames as a
# wide chunk file (rows=IDs, cols=new timepoints; no re-read), then trim the in-memory
# history to the last 2 frames (all hook calcs need only [-1]/[-2]). At the end the
# chunks are merged into the SAME counts_1.csv / SA_1.csv / Flux_1.csv format.
#
# Frame bookkeeping (hookInterval seconds per frame):
#   counts, SA : one frame per time 0,1,...,T  -> T+1 frames
#   flux       : one frame per time 1,...,T     -> T   frames (no t=0)
# =====================================================================================

def initCSVState(sim_properties):
    """Reset the streaming-output bookkeeping (call once before the CME run)."""
    sim_properties['csv_state'] = {'flush_idx': 0,
                                   'counts_written': 0,   # frames already flushed
                                   'sa_written': 0,
                                   'flux_written': 0}
    return None


def _chunk_dir(sim_properties):
    d = os.path.join(os.path.dirname(sim_properties['path']['counts']), '_csv_chunks')
    os.makedirs(d, exist_ok=True)
    return d


def _chunk_path(sim_properties, stream, idx):
    stem = os.path.basename(sim_properties['path'][stream])[:-4]   # strip '.csv'
    return os.path.join(_chunk_dir(sim_properties), '{0}.part{1:05d}.csv'.format(stem, idx))


def _write_wide_chunk(path, rowIDs, columns_times, value_lists):
    """value_lists[k] is the column of values for time columns_times[k] (len == len(rowIDs))."""
    df = pd.DataFrame()
    df['Time'] = rowIDs
    for t, vals in zip(columns_times, value_lists):
        df[t] = vals
    df.to_csv(path, index=False)


def flushChunk(sim_properties):
    """Write the frames accumulated since the last flush as wide chunk files
    (counts, SA, flux). No re-read of the main CSV. Call when the current time is
    a multiple of the output interval (and once more at the end of the run)."""

    st = sim_properties['csv_state']
    hookInterval = sim_properties['hookInterval']
    T = sim_properties['time_second'][-1]
    idx = st['flush_idx']

    # ---- counts (frames for times 0..T) ----
    counts = sim_properties['counts']
    n_total = int(round(T / hookInterval)) + 1
    num_new = n_total - st['counts_written']
    if num_new > 0:
        times = [(st['counts_written'] + i) * hookInterval for i in range(num_new)]
        cols = []
        for k in range(num_new):
            col = []
            for c in counts.values():
                # Some counts are not updated every second; fall back to the last value.
                try:
                    col.append(c[-num_new + k])
                except IndexError:
                    col.append(c[-1])
            cols.append(col)
        _write_wide_chunk(_chunk_path(sim_properties, 'counts', idx),
                          list(counts.keys()), times, cols)
        st['counts_written'] = n_total

    # ---- SA + volume_L (frames for times 0..T) ----
    SA = sim_properties['SA']
    vol = sim_properties['volume_L']
    num_new_sa = n_total - st['sa_written']
    if num_new_sa > 0:
        times = [(st['sa_written'] + i) * hookInterval for i in range(num_new_sa)]
        cols = []
        for k in range(num_new_sa):
            col = [s[-num_new_sa + k] for s in SA.values()]
            col.append(vol[-num_new_sa + k])
            cols.append(col)
        _write_wide_chunk(_chunk_path(sim_properties, 'SA', idx),
                          list(SA.keys()) + ['volume_L'], times, cols)
        st['sa_written'] = n_total

    # ---- flux (frames for times 1..T; offset by one, no t=0) ----
    fluxes = sim_properties.get('fluxes')
    if fluxes:
        n_total_flux = int(round(T / hookInterval))
        num_new_flux = n_total_flux - st['flux_written']
        if num_new_flux > 0:
            times = [(st['flux_written'] + 1 + i) * hookInterval for i in range(num_new_flux)]
            cols = []
            for k in range(num_new_flux):
                col = [f[-num_new_flux + k] for f in fluxes.values()]
                cols.append(col)
            _write_wide_chunk(_chunk_path(sim_properties, 'flux', idx),
                              list(fluxes.keys()), times, cols)
            st['flux_written'] = n_total_flux

    st['flush_idx'] += 1
    return None


def trimHistory(sim_properties, keep=2):
    """Trim the per-hook histories to the last `keep` frames (default 2: hook cost
    calcs need [-1] and [-2]). time_second is the time axis and is NOT trimmed."""
    for c in sim_properties['counts'].values():
        del c[:-keep]
    for s in sim_properties['SA'].values():
        del s[:-keep]
    del sim_properties['volume_L'][:-keep]
    fluxes = sim_properties.get('fluxes')
    if fluxes:
        for f in fluxes.values():
            del f[:-keep]
    return None


def mergeChunks(sim_properties):
    """Merge the per-flush wide chunks into the final counts_1.csv / SA_1.csv /
    Flux_1.csv (identical format to the legacy writers), then delete the chunks."""
    d = _chunk_dir(sim_properties)
    for stream in ['counts', 'SA', 'flux']:
        stem = os.path.basename(sim_properties['path'][stream])[:-4]
        parts = sorted(glob.glob(os.path.join(d, stem + '.part*.csv')))
        if not parts:
            continue
        merged = pd.read_csv(parts[0])
        for p in parts[1:]:
            merged = pd.concat([merged, pd.read_csv(p).drop(columns=['Time'])], axis=1)
        merged.to_csv(sim_properties['path'][stream], index=False)
        for p in parts:
            os.remove(p)
    try:
        os.rmdir(d)
    except OSError:
        pass
    return None


def writeCountstoCSV(restartNum, sim_properties):
    """

    Description:  Write the trajectories of all species into CSV file per restartCME
                Also filling the possible shorten trajectories in sim_properties['counts'] due to the jumping of large waiting time by the last value of each trajectory 
    """
    
    countsfile_path = sim_properties['path']['counts']

    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    currenttime_second = sim_properties['time_second'][-1]


    if restartNum == 0:
        particleDF = pd.DataFrame()

        spec_IDs = []

        for specID in sim_properties['counts'].keys():
            spec_IDs.append(specID)
        
        particleDF['Time'] = spec_IDs

        hookMoments = np.arange(restartNum*restartInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        for hookMoment in hookMoments:
            
            # location of trajectory in countsDic
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:
                print('Time moment {0} second is not covered between {1} second and {2} second due to a long reaction waiting time'.
                      format(hookMoment, restartNum*restartInterval,(restartNum+1)*restartInterval) )
                for count in sim_properties['counts'].values():
                    count.append(count[-1])
                print('The counts of species at time {0} second are repeated by count at time {1} second '.format(hookMoment, currenttime_second))

            new_counts = []
            for species, count in sim_properties['counts'].items():
                # 1) Certain counts are not updated during the simulation so use try-expect to avoid list index out of range
                try:
                    new_counts.append(count[i_loc])
                except:
                    new_counts.append(count[-1])

            particleDF[hookMoment] = new_counts
        particleDF.to_csv(countsfile_path,index=False)

    else:
        hookMoments = np.arange(restartNum*restartInterval+hookInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        particleDF = pd.read_csv(countsfile_path)

        for hookMoment in hookMoments:
    
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:
                print('Time moment {0} second is not covered between {1} second and {2} second due to a long reaction waiting time'.
                      format(hookMoment, restartNum*restartInterval,(restartNum+1)*restartInterval) )

                for count in sim_properties['counts'].values():
                    count.append(count[-1])
                print('The counts of species at time {0} second are repeated by count at time {1} second '.format(hookMoment, currenttime_second))

            new_counts = []
            for count in sim_properties['counts'].values():
                # Certain counts are not updated per second so use try-expect to avoid list index out of range
                try:
                    new_counts.append(count[i_loc])
                except:
                    new_counts.append(count[-1])

                
            particleDF[hookMoment] = new_counts

        particleDF.to_csv(countsfile_path,index=False)

    return None

def writeSAtoCSV(restartNum, sim_properties):
    """
    
    Description: Write the trajectories of Surface area and volume into CSV file per restartCME
        Also filling the possible shorten trajectories in sim_properties['SA'] due to the jumping of large waiting time by the last value of each trajectory 
    """

    SAfile_path = sim_properties['path']['SA']


    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    currenttime_second = sim_properties['time_second'][-1]

    if restartNum == 0:
        SADF = pd.DataFrame()

        IDs = []

        for specID in sim_properties['SA'].keys():
            IDs.append(specID)
        
        IDs.append('volume_L')

        SADF['Time'] = IDs

        hookMoments = np.arange(restartNum*restartInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        for hookMoment in hookMoments:
            
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:

                # print('Time moment {0} second is not covered between {1} second and {2} second due to a long reaction waiting time'.
                      # format(hookMoment, restartNum*restartInterval,(restartNum+1)*restartInterval) )
                #sim_properties['time_second'].append(hookMoment)
                #print('{0} second is appened to the hookMoment')
                for number in sim_properties['SA'].values():
                    number.append(number[-1])
                sim_properties['volume_L'].append(sim_properties['volume_L'][-1])

                print('The SA and volume at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))
            # 0, 1, 2, ..., 60 second
            new_numbers = []
            for numbers in sim_properties['SA'].values():
                new_numbers.append(numbers[i_loc])

            new_numbers.append(sim_properties['volume_L'][i_loc])

            SADF[hookMoment] = new_numbers

        SADF.to_csv(SAfile_path,index=False)

    else:
        particleDF = pd.read_csv(SAfile_path)

        hookMoments = np.arange(restartNum*restartInterval+hookInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

        for hookMoment in hookMoments:
            # 61, 62, ..., 120 second, ...
            
            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:

                for number in sim_properties['SA'].values():
                    number.append(number[-1])
                sim_properties['volume_L'].append(sim_properties['volume_L'][-1])

                print('The SA and volume at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))
            
            new_numbers = []
            for numbers in sim_properties['SA'].values():
                new_numbers.append(numbers[i_loc])

            new_numbers.append(sim_properties['volume_L'][i_loc])

            particleDF[hookMoment] = new_numbers

        particleDF.to_csv(SAfile_path,index=False)


    return None



def round_sig(x, sig=2):
    negative = False
    if x < 0:
        negative = True
    x = abs(x)
    if negative:
        return -1*round(x, sig-int(floor(log10(abs(x))))-1)
    elif x==0.0:
        return 0.0
    else:
        return round(x, sig-int(floor(log10(abs(x))))-1)


def writeFluxtoCSV(restartNum, sim_properties):
    """
    Called after the finish of one CME simulation

    Description: Write the fluxes through ODE reactions per restartCME
    # Since the kinetic parameters remain during the simulation, we can generate the fluxes based on one odecell object
    # The difficuly for per CMEInterval is how to get the concentration list of metabolites
    # We will create a new subdictionary called sim_properties['conc'] to store the concetrations list
    # the conc subdictionary starts from hookInterval seconds not 0 seconds
    """
    # Initialize odecell again
    # odemodel = ODE.initModel(sim_properties)

    fluxesDict = sim_properties['fluxes']

    Fluxfile_path = sim_properties['path']['flux']

    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    currenttime_second = sim_properties['time_second'][-1]

    hookMoments = np.arange(restartNum*restartInterval+hookInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

    if restartNum == 0:
        FluxDF = pd.DataFrame()

        rxnIDs = fluxesDict.keys()

        FluxDF['Time'] = rxnIDs

        for hookMoment in hookMoments:
            # 1, 2, 3, ... restartInterval
            
            i_loc = int(hookMoment/hookInterval)
            if hookMoment > currenttime_second:
                for flux in fluxesDict.values():

                    flux.append(flux[-1])
                print('The fluxes at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))
            
            flux_hookMoment = [flux[i_loc-1] for flux in fluxesDict.values()]

            FluxDF[hookMoment] = flux_hookMoment

        FluxDF.to_csv(Fluxfile_path,index=False)

    else:
        FluxDF = pd.read_csv(Fluxfile_path)

        for hookMoment in hookMoments:

            i_loc = int(hookMoment/hookInterval)

            if hookMoment > currenttime_second:
                for flux in fluxesDict.values():

                    flux.append(flux[-1])
                print('The fluxes at time {0} second are repeated by that at time {1} second'.format(hookMoment, currenttime_second))

            flux_hookMoment = [flux[i_loc-1] for flux in fluxesDict.values()]

            FluxDF[hookMoment] = flux_hookMoment

        FluxDF.to_csv(Fluxfile_path,index=False)

    return None


def appendHookMoments(restartNum, sim_properties):
    """
    Description:
        To make time moments in sim_properties complete by appending time points
        solve the issue that the possible jumping of CME time over several hookIntervals 
    """

    restartInterval = sim_properties['restartInterval']
    
    hookInterval = sim_properties['hookInterval']

    hookMoments = np.arange(restartNum*restartInterval,(restartNum+1)*restartInterval+hookInterval, hookInterval)

    currenttime_second = sim_properties['time_second'][-1]

    for hookMoment in hookMoments:

        if hookMoment > currenttime_second:

            sim_properties['time_second'].append(hookMoment)
            
            print('{0} second is appened to the simulationTime'.format(hookMoment))

    return None