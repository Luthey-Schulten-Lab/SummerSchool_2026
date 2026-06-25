"""
Author: Enguang Fu
Date: March 2024

Main Scirpt for CMEODE Hybrid Whole Cell Simulation
"""

###################################
#####   Import ####
###################################

# pyLM to construct the CME system
from pyLM import *

# 
import numpy as np
import os, sys
import time as phys_time
from Bio import SeqIO

from datetime import datetime
import argparse

from mpi4py import MPI

# User defined Python scripts
import initiation as IC
import rxns_CME
import communicate
import species_counts
import hookSolver_CMEODE
import filesaving
import hook_CMEODE



###################################
#####   Parsing Parameters ####
###################################


comm = MPI.COMM_WORLD
rank = comm.Get_rank()+1
size = comm.Get_size()


# Arguments needed: input_dir, simulationType, simulation time length, restartInterval, hookInterval, writeInterval, outputfolder
ap = argparse.ArgumentParser()


ap.add_argument('-in','--input_dir',required=True)
ap.add_argument('-st','--simType',required=True)
ap.add_argument('-t','--simTime',required=True)
# -o/--outputInterval is the CSV output cadence (the CME no longer restarts);
# -rs/--restartInterval kept as a backward-compatible alias.
ap.add_argument('-o','--outputInterval','-rs','--restartInterval',required=True, dest='restartInterval')
ap.add_argument('-hi', '--hookInterval', required = True)
ap.add_argument('-f', '--output_dir', required= True)

args = ap.parse_args()

input_dir = args.input_dir

output_dir = args.output_dir

# Check if the directory exists
if not os.path.exists(output_dir):
        # The directory does not exist, create it
        os.makedirs(output_dir)


logfile = 'log_{0}.txt'.format(rank)

logfile_path = os.path.join(output_dir, logfile)

# Redirect the output and possible error to a log file
log = open(logfile_path,'w')

original_stdout = sys.stdout
original_stderr = sys.stderr

sys.stdout = log
sys.stderr = log


print('The PID is ' + str(os.getpid()))

try:
    
    simTime = float(args.simTime)
    restartInterval = float(args.restartInterval)
    hookInterval = float(args.hookInterval)

except:
    
    sys.exit("Error: Please enter simulation time (-t), restart time (-rs), hook interval (-hi), and write interval (-wi) as the valued time in seconds.")

if simTime % restartInterval == 0:
        if restartInterval % hookInterval == 0:
                    if (str(args.simType) == "cme-ode"):
                           None
                    else:
                                sys.exit("Error: Enter 'cme-ode' for a hybrid CME-ODE simulation")
        else:
                sys.exit('Error: Please enter restart time (-rs) as integer multiples of hook interval (-hi)')
else:
        sys.exit("Error: Please enter simulation time (-t) as integer multiples of restart time (-rs)")


# -rs is now the CSV OUTPUT interval (NO CME restart): the CME runs ONCE for the
# whole simulation time and the per-hook histories are flushed to CSV (and trimmed
# to the last 2 frames) every outputInterval seconds, from inside the hook.
# update_rateconstants already refreshes CME rates every hook, so the periodic
# CME restart that this used to drive is no longer needed.
outputInterval = restartInterval
outputNums = int(simTime/outputInterval)

print('The simulation time is {0} seconds; CSV output every {1} seconds (single CME run, no restart).'.format(simTime, outputInterval))
print('The hook interval is {0} seconds and write interval to the LM file {1} seconds'.format(hookInterval, hookInterval))
print('CME runs once; CSV will be flushed {0} times'.format(outputNums))


print('*******************************************************************************')


######################################
#####   Initialize sim_properties ####
######################################

start_time = datetime.now()

print('The simulation starts at '+ str(start_time))

# Initialize sim_properties, a dictionary record the genome information, trajectories, and other constants

sim_properties = {}

sim_properties['input_dir'] = input_dir
print(f"The input directory is {input_dir}")

# Input Gene bank file syn3A.gb
genomeFile3A =  sim_properties['input_dir'] + 'syn3A.gb'
genome3A = next(SeqIO.parse(genomeFile3A, "gb"))

# Convert genebank file syn3A.gb into a multilayer dictionanry genome
sim_properties['genome'], sim_properties['genome_length'] = IC.mapDNA(genome3A)

# Initialize the constants
IC.initializeConstants(sim_properties)

# Output the corrected counts for proteins in complexes into a Excel sheet
IC.outPtnCorrInitCounts(sim_properties, rank)

# Map the locusNum to gene and intergen region
IC.getlocusNumtoGeneSeq(sim_properties, genome3A)

# Initialize the counts and concentrations of metabolites
IC.initializeMetabolitesCounts(sim_properties)
IC.initializeMediumConcs(sim_properties)
IC.initializeProteinMetabolitesCounts(sim_properties)

# Initialize the traces and counts of cost booking species
IC.initializeCosts(sim_properties)

# Initialize to construct the reaction maps
IC.getReactionMap(sim_properties)

# The names of three csv files e.g. counts_1.csv, SA_1.csv, Flux_1.csv
countsCSV = 'counts_{0}.csv'.format(rank)
countsCSV_path = os.path.join(output_dir,countsCSV)

SACSV = 'SA_{0}.csv'.format(rank)
SACSV_path = os.path.join(output_dir, SACSV)

fluxCSV = 'Flux_{0}.csv'.format(rank)
fluxCSV_path = os.path.join(output_dir, fluxCSV)

sim_properties['path'] = {'counts': countsCSV_path, 'SA': SACSV_path, 'flux':fluxCSV_path}

# sim_properties['time_second] will append the hook moments until hits the simulation time length
sim_properties['time_second'] = [int(0)]

# simulaion parameters
sim_properties['restartInterval'] = restartInterval

sim_properties['hookInterval'] = hookInterval

sim_properties['restartNum'] = []

# Single-CME-run output settings: full simulation time + CSV flush cadence.
sim_properties['simTime'] = simTime
sim_properties['outputInterval'] = outputInterval

# Initialize the append-only streaming-CSV bookkeeping.
filesaving.initCSVState(sim_properties)

# PROFILING: wall-time accumulators (seconds) for the restart-cost breakdown.
sim_properties['prof'] = {'cme_setup': 0.0, 'cme_save': 0.0, 'cme_runsolver': 0.0,
                          'ode_init': 0.0, 'ode_sim': 0.0,
                          'hook_total': 0.0, 'hook_count': 0}


###################################
#####  Perform the Simulation  ####
###################################


# NO CME restart: build and run the CME exactly ONCE for the whole simulation
# time. The loop is kept (single iteration) to minimize structural change; CSV
# output now happens periodically from inside the hook (see hookSolver_CMEODE).
for restartNum in range(0, 1):

    sim_properties['restartNum'].append(restartNum)
    print('####################################################################################')

    print('Start the single CME run over 0 to {0} seconds (CSV output every {1} s)'.format(simTime, outputInterval))

    # Create one instance of Class CMESimulation
    sim=CME.CMESimulation(name="CMEODE_" + str(restartNum+1))

    # Set the parameters. writeInterval, hookInterval and simulationTime of the CME
    sim = IC.initializeCME(sim, restartNum, sim_properties)

    # Run the WHOLE simulation in one CME (override the per-restart simulationTime).
    sim.setSimulationTime(simTime)

    # define the species, reactions and particles numbers for every newly start simulation
    ini_start = phys_time.time()

    rxns_CME.addGeneticInformationReactions(sim, sim_properties, genome3A)

    IC.addGeneticInformationSpeciesCounts(sim,sim_properties)

    sim_properties['prof']['cme_setup'] += phys_time.time() - ini_start
    print(f'Initiation of CME reactions and counts in {phys_time.time() - ini_start:.3f} seconds')

    # Initialize CME_count_array after setting up the entire CME system
    CME_count_array = species_counts.SpeciesCounts(sim)

    hook = hookSolver_CMEODE.MyOwnSolver()
    
    hook.initializeSolver(sim, CME_count_array, sim_properties, genome3A)

    # Initialize the membrane based on lipid and protein counts
    if restartNum == 0:
        IC.initializeMembrane(sim_properties)
    
    # LMfilename of the current simulation
    LMfilename = os.path.join(output_dir, 'CME_ODE_{0}_{1}.lm'.format(rank, restartNum))

    # Remove lm file from previous simulation and Save new CME simulation LM HDF5 file
    os.system("rm -rf %s"%(LMfilename))
    save_start = phys_time.time()
    sim.save(LMfilename)
    sim_properties['prof']['cme_save'] += phys_time.time() - save_start
    print(f'LM file Saved in {phys_time.time() - save_start:.3f} seconds')

    print(sim_properties['locusNumtoGeneSeq'])
    
    if restartNum == 0:

        # print(sim.initial_counts)

        print('Genetic Information Reactions Numbers in different subsystem:')
        print(sim_properties['rxns_numbers'])

        print("Total reactions in sim_properties['rxns_numbers'] are {0}".format(sum(number_rxns[0] for number_rxns in sim_properties['rxns_numbers'].values()) ))
        
        print('rxns_map')
        print(f"{sim_properties['rxns_map']}")

        # Print CME system information
        print('{0} Species in CME simulation object sim.particleMap'.format(len(sim.particleMap)))

        print(sim.particleMap)
        print('self.initial_counts')
        print(f'{sim.initial_counts}')
        print('{0} Reactions in CME simulation object sim.reactions'.format(len(sim.reactions)))
        
        for rxn in sim.reactions:
               print(rxn) # For debugging

        # print(sim_properties['rxns_map'])
        print(sim_properties['counts'])


    runCMEODE_start = phys_time.time()

    # run CME simulation with hook of one single replicate
    runsolver_start = phys_time.time()
    sim.runSolver(filename=LMfilename, solver=hook, replicates=1)
    sim_properties['prof']['cme_runsolver'] += phys_time.time() - runsolver_start

    # Now the CME simulation ends, we use updateCMEcountsFile to read in the counts and pass them to ODE to do the last second's ODE simulation

    print(f'HookSimulation is called at {simTime:.3f} second')

    communicate.updateCMEcountsFile(sim, sim_properties, LMfilename)

    # run hook again to make sure we run N instead of N-1 times hook (final second)
    hook_CMEODE.hook_CMEODE(sim_properties, genome3A)

    CSVfilesaving_start = phys_time.time()

    # Final flush of any frames since the last periodic (in-hook) flush, then trim.
    filesaving.flushChunk(sim_properties)
    filesaving.trimHistory(sim_properties, keep=2)

    print(f'Final CSV chunk flushed in {phys_time.time() - CSVfilesaving_start:.3f} seconds')


    # Remove the just finished lm file
    os.system("rm -rf %s"%(LMfilename))

    print(f'Finish the single CME run (0 to {simTime} seconds) in {phys_time.time() - runCMEODE_start:.3f} seconds')

    print('####################################################################################')

# Merge the per-flush CSV chunks into the final counts/SA/Flux CSVs (same format).
merge_start = phys_time.time()
filesaving.mergeChunks(sim_properties)
print(f'CSV chunks merged into final counts/SA/Flux files in {phys_time.time() - merge_start:.3f} seconds')



###################################
#####   Close the Simulation ####
###################################

end_time = datetime.now()

print('The simulation ends at {0}'.format(end_time))

print('Time (hour:minutes:seconds) taken to finish the simulation: {0}'.format(end_time - start_time))

# PROFILING SUMMARY: restart-cost breakdown over the whole run.
_p = sim_properties['prof']
_cme_pure = _p['cme_runsolver'] - _p['hook_total']          # CME C++ solver, hooks removed
_hook_other = _p['hook_total'] - _p['ode_init'] - _p['ode_sim']  # communicate/rate-update overhead
_wall = (end_time - start_time).total_seconds()
print('==================== PROFILING SUMMARY ====================')
print(f'CME runs=1 (no restart)  CSV flushes={outputNums}  hooks(in-solver)={_p["hook_count"]}  wall={_wall:.1f}s')
print(f'  CME setup (addGeneticInfo+counts) : {_p["cme_setup"]:8.2f} s   [ONCE]')
print(f'  CME save  (sim.save .lm)          : {_p["cme_save"]:8.2f} s   [ONCE]')
print(f'  CME simulation (pure C++ solver)  : {_cme_pure:8.2f} s')
print(f'  ODE initialize (initModel)        : {_p["ode_init"]:8.2f} s')
print(f'  ODE simulation (makeSolver+runODE): {_p["ode_sim"]:8.2f} s')
print(f'  hook other (communicate/rates)    : {_hook_other:8.2f} s')
print(f'  --- runSolver wall (CME+hooks)    : {_p["cme_runsolver"]:8.2f} s ; hook_total {_p["hook_total"]:.2f} s')
print('===========================================================')

sys.stdout = original_stdout
sys.stderr = original_stderr

# print out in the terminal
print('The simulation ends at {0}'.format(end_time))

print('Time taken to finish the simulation: {0}'.format(end_time - start_time))

