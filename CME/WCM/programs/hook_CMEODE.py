"""
Author: Enguang Fu

Date: March 2024

the operations when do hooking/communicating between CME/ODE
"""

import time as phys_time
import numpy as np

import communicate
import integrate
import rxns_ODE as ODE
import filesaving
import initiation as IC


def hook_CMEODE(sim_properties, genome3A):
    """
    Input: sim_properties, genome3A
    
    Description: 1) calculate the cost of nucleotides and pass the extracted counts of ODE species to odecell;
                 2) Do ODE simulation
                 3) Record the ODE counts into counts and conc dictionaries in sim_properties and update ODE result to CME
                 3) update the Surface area
    """
    
    # save the transient end flux or not
    boolean_end = False
    
    # calculate the cost of gene expression
    communicate.calculateCosts(sim_properties, genome3A, sim_properties['locusNumtoIndex'])

    # Communicate the cost of monomers; Record Shortage
    communicate.communicateCostsToMetabolism(sim_properties)

    # Check whether negative count before ODE run
    communicate.checkbeforeODE(sim_properties)

    iniode_start = phys_time.time()

    # Initialize/construct the ODEs using odecell
    # The counts in the counts dictionary are passed to ODE solver as initial values
    odemodel = ODE.initModel(sim_properties)

    # Create flux dictionary to store the reaction fluxes
    if sim_properties['time_second'][-1] == 0:
        IC.initializeFluxes(sim_properties, odemodel, boolean_end=boolean_end)
    
    print(f'ODE object Initialized in {phys_time.time() - iniode_start:.3f} seconds')

    runode_start = phys_time.time()

    # Convert odemodel to right hand side of the differential equations
    solver = integrate.noCythonSetSolver(odemodel)

    odelength = sim_properties['hookInterval']

    concs, flux_avg, flux_end = integrate.runODE(solver, odemodel, odelength)

    print(f'ODE simulation of {odelength} second Finished in {phys_time.time()- runode_start:.3f} seconds')
    
    if phys_time.time()- runode_start > 10:
        print(f"WARNING: Simulation Time {sim_properties['time_second'][-1]} ODE Integral takes {phys_time.time()- runode_start:.3f} seconds")
        
    # save the flues into sim_properties['fluxes']
    communicate.updateFluxes(sim_properties, odemodel, flux_avg, flux_end, boolean_end=boolean_end)
    
    # Save the concentrations directly from ODE to conc dictionary 
    # communicate.saveConc(sim_properties, odeResults, odemodel)
    
    # update the counts of ODE species into counts dictionary
    communicate.updateODEcounts(sim_properties, concs, odemodel)

    # update the Surface Area and volume
    communicate.updateSA(sim_properties)
    
    # forward time by hookInterval
    sim_properties['time_second'].append(sim_properties['time_second'][-1] + sim_properties['hookInterval'])

    checkuint_32(sim_properties)

    return None


def checkuint_32(sim_properties):
    """

    Description: Check the data type of counts stored in sim_properties per communication step
    """
    countsDic = sim_properties['counts']

    for specie, count in countsDic.items():

        value = count[-1]
        if not (value >= 0 and (type(value) == int or type(value) == np.uint32)):
            print(f"Time: {sim_properties['time_second'][-1]} {specie} {value} {type(value)}")

    return None