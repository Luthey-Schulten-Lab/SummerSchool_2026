"""
Authors: Zane Thornburg, Enguang Fu

A file to define all of the reactions in the Metabolism
"""


import os
import sys

import odecell

import libsbml
import pandas as pd
import numpy as np

import integrate

from collections import defaultdict

#########################################################################################
def initModel(sim_properties):
    """
    Initiate the model object and pass in the CME species counts

    Arguments:

    sim_properties['counts'] (particle map): The CME particle Map

    Returns:

    upModel (odecell model object): The updated ODE kinetic model for simulation

    Acceleration (Phase 1, build-once): the ODE model topology (reactions, rate
    forms, kcat, KM, stoichiometry, medium params, metabolite set) is STATIC for
    the whole run; only metabolite concentrations and enzyme concentrations change
    per hook. So the model is built ONCE (reading Excel/SBML once) and cached in
    sim_properties['_ode_cache']; subsequent hooks only refresh the dynamic values
    on the cached model. This removes the per-hook file-I/O + model construction
    (~0.65 s/hook) with no change to the numerical result.

    Phase 2 (Cython): on the first hook the RHS is additionally compiled ONCE
    into a Cython functor (enzyme concentrations routed through self.params[] so
    they can be updated each hook with no recompile -- see makeSolver). The
    integration then runs the compiled C RHS (~5x faster). Controlled by
    sim_properties['ode_use_cython'] (default True); on any build failure it
    latches to the no-Cython cached path (Fallback 1) for the rest of the run.

    Tiers:
      - sim_properties['ode_no_cache'] = True  -> Fallback 2: original per-hook
        rebuild (bit-for-bit reproduction / debugging).
      - sim_properties['ode_use_cython'] = False -> Phase 1 only (cached model,
        no-Cython codegen each hook).
      - default -> Phase 1 cache + Phase 2 Cython functor.
    """

    # Fallback 2: original behaviour, full rebuild every call.
    if sim_properties.get('ode_no_cache', False):
        return _buildModel(sim_properties, cache=None)

    cache = sim_properties.get('_ode_cache')

    if cache is not None:
        # Subsequent hooks: reuse the cached model, refresh dynamic values only.
        # Metabolite y0 is needed by BOTH paths (it is the integrator state).
        _refreshMetaboliteInitVals(sim_properties, cache)
        # Enzyme concs are baked into the RHS only on the no-Cython path; the
        # Cython path injects them via the param vector in makeSolver instead.
        if not cache.get('cython', False):
            _refreshEnzymeParams(sim_properties, cache)
        return cache['model']

    # First hook: build once and record the metadata needed to refresh later.
    cache = {'enzyme_specs': [], 'base_forms': {}, 'enzyme_spec_by_rxnID': {}}
    model = _buildModel(sim_properties, cache=cache)
    cache['model'] = model
    sim_properties['_ode_cache'] = cache

    # Phase 2: try to compile the Cython functor once. On failure, latch to the
    # no-Cython cached path (the model is rebuilt clean so enzyme baking works).
    cache['cython'] = False
    if sim_properties.get('ode_use_cython', True):
        try:
            _setupCythonFunctor(sim_properties, cache)
            cache['cython'] = True
            print("ODE: Cython functor compiled and cached (build-once).")
        except Exception as exc:
            print(f"ODE: Cython functor setup FAILED ({exc!r}); "
                  f"falling back to no-Cython cached path for the rest of the run.")
            cache['model'] = _buildModel(sim_properties, cache=None)  # clean rebuild
            cache['cython'] = False

    return cache['model']
#########################################################################################


#########################################################################################
def makeSolver(sim_properties, odemodel):
    """Return the ODE RHS solver for this hook.

    Cython path: re-instantiate the cached compiled functor with a fresh param
    vector (onoff=1, enzyme concs recomputed from current counts) -- a cheap
    array copy, no recompile. No-Cython path: the cheap per-hook codegen
    (integrate.noCythonSetSolver), which re-bakes the enzyme params updated by
    _refreshEnzymeParams.
    """

    cache = sim_properties.get('_ode_cache')

    if cache is not None and cache.get('cython', False):
        vec = _buildParamVector(sim_properties, cache)
        return cache['functor_class'](np.asarray(vec, dtype=np.double))

    return integrate.noCythonSetSolver(odemodel)
#########################################################################################


#########################################################################################
def _buildModel(sim_properties, cache=None):
    """Construct the odecell model from scratch (reads Excel/SBML).

    When `cache` is provided, refresh metadata (enzyme specs, protein base-form
    recipes) is recorded into it so the model can be refreshed in place later.
    """

    # Initialize the ODECell model
    model = odecell.modelbuilder.MetabolicModel()

    zeroOrderOnOff = '$onoff * $K'

    model.zeroOrderOnOff = odecell.modelbuilder.RateForm(zeroOrderOnOff)

    model.updateAvailableForms()

    # Set verbosity outputs to zero for now to improve performance
    model.setVerbosity(0)

    # Define Rxns and pass in the Particle Map containing enzyme concentrations
    model = defineRxns(model, sim_properties, cache=cache)

    return model
#########################################################################################


#########################################################################################
def _setupCythonFunctor(sim_properties, cache):
    """Compile the Cython functor ONCE and cache the class + param-vector layout.

    'Enzyme' params are promoted to optimizable (lb != ub) so stock odecell's
    prepareFunctor routes them through self.params[]; the functor is then
    instantiated each hook with updated enzyme values (no recompile). The build
    runs in a per-process directory to avoid the MPI cwd collision (odecell
    writes cythonCompiledFunctions.pyx / setup_tmp.py into the cwd).
    """

    model = cache['model']

    # Promote every 'Enzyme' reaction parameter to optimizable so it routes
    # through self.params[] in the compiled functor.
    _, _, rxnParam = model.getParameters()
    for rxnID, parDict in rxnParam.items():
        if 'Enzyme' in parDict:
            parDict['Enzyme'].ub = 1.0e6

    # Capture the opt-space layout BEFORE prepareFunctor replaces the opt-param
    # values with 'self.params[i]' strings. currParVals is the correctly-ordered,
    # correctly-sized numeric vector (this is what 4DWCM got wrong).
    optSpace, currParVals = model.getOptSpace()
    fields = set(o.field for o in optSpace)
    if not fields <= {'onoff', 'Enzyme'}:
        raise RuntimeError(f"unexpected optimizable params in ODE model: {fields - {'onoff', 'Enzyme'}}")

    cache['param_template'] = [float(v) for v in currParVals]
    # slots that must be refreshed each hook: (vector_index, rxnID)
    cache['enzyme_slots'] = [(i, optSpace[i].indx) for i in range(len(optSpace))
                             if optSpace[i].field == 'Enzyme']

    solver = odecell.solver.ModelSolver(model)
    solver.prepareFunctor()

    # Silence the benign NumPy "deprecated API" #warning (category -Wcpp) emitted
    # by gcc when compiling the Cython-generated C. The build runs in a subprocess
    # spawned by buildCall, which inherits os.environ, so appending -Wno-cpp to
    # CFLAGS (preserving conda's existing flags) keeps the build log clean. Using
    # -Wno-cpp rather than -DNPY_NO_DEPRECATED_API avoids any risk of build errors
    # if the generated C touches the deprecated API.
    if '-Wno-cpp' not in os.environ.get('CFLAGS', ''):
        os.environ['CFLAGS'] = (os.environ.get('CFLAGS', '') + ' -Wno-cpp').strip()

    # Build in a per-process dir so concurrent MPI ranks do not clobber each
    # other's cythonCompiledFunctions.{pyx,c,so} / setup_tmp.py. odecell writes
    # those files into the cwd and then imports the module, so we both chdir into
    # the build dir AND put it on sys.path (the script dir, not '', is sys.path[0]).
    out_dir = os.path.dirname(sim_properties.get('path', {}).get('counts', '')) or '.'
    build_dir = os.path.abspath(os.path.join(out_dir, '_ode_cython_build_{0}'.format(os.getpid())))
    os.makedirs(build_dir, exist_ok=True)
    cwd = os.getcwd()
    sys.path.insert(0, build_dir)
    try:
        os.chdir(build_dir)
        solver.buildCall(odeint=False, useJac=False, cythonBuild=True,
                         functor=True, verbose=0)
    finally:
        os.chdir(cwd)
        try:
            sys.path.remove(build_dir)
        except ValueError:
            pass

    cache['functor_class'] = solver.functor
    return None
#########################################################################################


#########################################################################################
def _buildParamVector(sim_properties, cache):
    """Build the functor parameter vector for this hook: start from the build-time
    template (keeps onoff = 1) and overwrite the enzyme slots with concentrations
    recomputed from the current counts (same numbers the no-Cython refresh uses)."""

    vec = list(cache['param_template'])
    spec_by_rxnID = cache['enzyme_spec_by_rxnID']
    for slot, rxnID in cache['enzyme_slots']:
        EnzymeStr, GPRrule = spec_by_rxnID[rxnID]
        vec[slot] = _enzymeConcFromSpec(EnzymeStr, GPRrule, sim_properties)
    return vec
#########################################################################################


#########################################################################################
def _refreshEnzymeParams(sim_properties, cache):
    """No-Cython path: update enzyme concentrations on the cached model (re-baked
    by the per-hook noCythonSetSolver codegen)."""

    model = cache['model']
    for rxnIndx, EnzymeStr, GPRrule in cache['enzyme_specs']:
        model.addParameter(rxnIndx, "Enzyme",
                           _enzymeConcFromSpec(EnzymeStr, GPRrule, sim_properties))
    return None
#########################################################################################


#########################################################################################
def _refreshMetaboliteInitVals(sim_properties, cache):
    """Update metabolite initial concentrations (the integrator state vector y0)
    on the cached model via setInitVal, so model.getInitVals() returns fresh
    values. Needed by both the Cython and no-Cython paths."""

    model = cache['model']
    counts = sim_properties['counts']
    base_forms = cache['base_forms']
    for metID, idx in model.getMetDict().items():
        if metID in base_forms:
            ptnID, otherForms = base_forms[metID]
            baseFormCount = int(counts[ptnID][-1] - sum(counts[m][-1] for m in otherForms))
            conc = partTomM(baseFormCount, sim_properties)
        else:
            conc = partTomM(counts[metID][-1], sim_properties)
        model.getMetList(idx).setInitVal(conc)
    return None
#########################################################################################


#########################################################################################
def defineRxns(model, sim_properties, cache=None):

    model = addProteinMetabolites(model, sim_properties, cache=cache)

    model = defineRandomBindingRxns(model, sim_properties, cache=cache)

    model = defineNonRandomBindingRxns(model, sim_properties)
    
    # Other Random Bindig reactions are GLCK and GLCT reactions that converts extracellular glucose into g6p.
    # model = defineOtherRandomBindingReactions(model, sim_properties)
    
    return model
#########################################################################################


#########################################################################################
def partTomM(particles, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's

    conc = (particles*1000.0)/(NA*sim_properties['volume_L'][-1])

    return conc
#########################################################################################


#########################################################################################
def reptModel(model):
    """
    Report on the constructed hybrid model - but probably would only want to do after the first time step

    Arguments: 
    model (model obj.): The ODE Kinetic Model

    Returns:

    None
    """

    dictTypes = defaultdict(int)
    typeList = ["Transcription","Translation","Degradation"]

    for rxn in model.getRxnList():
        
        if rxn.getResult():
            # If an explicit result has been set, this is a dependent reaction.
            dictTypes["Dependent reactions"] += 1
            continue
        
        for rxntype in typeList:
            if rxntype in rxn.getID():
                dictTypes[rxntype] += 1

                
    #print( "There are {} ODEs in the model:".format(len(model.getRxnList())) )

    outList = list(dictTypes.items())
    outList.sort(key=lambda x: x[1], reverse=True)
    for key,val in outList:
        print("{:>20} :   {}".format(key,val) )
        return 0
    
    return None
#########################################################################################

# Central, Nucleotide, Lipid, Cofactor and Transprot Reactions were defined.
#########################################################################################
def defineRandomBindingRxns(model, sim_properties, cache=None):
    """
    Define all of the reactions and rateforms needed for the current module to an existing module.

    """
    
    params_file = sim_properties['kinetic_params_path']
    
    central_params = pd.read_excel(params_file, sheet_name='Central')
    nucleotide_params = pd.read_excel(params_file, sheet_name='Nucleotide')
    lipid_params = pd.read_excel(params_file, sheet_name='Lipid')
    cofactor_params = pd.read_excel(params_file, sheet_name='Cofactor')
    transport_params = pd.read_excel(params_file, sheet_name='Transport')
    
    metabolism_params = pd.concat([central_params, nucleotide_params, lipid_params, cofactor_params, transport_params], ignore_index=True) #, transport_params
    
    reaction_list = []

    for row, item in metabolism_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])


    # The .xml file for flux balance analysis
            
    sbmlFile = "../input_data/Syn3A_updated.xml"

    docSBML = libsbml.readSBMLFromFile(sbmlFile)
    modelSBML = docSBML.getModel()

    speciesNames = [spc.name for spc in modelSBML.getListOfSpecies()]
    speciesNamesLower = [x.lower() for x in speciesNames]
    speciesIDs = [spc.id for spc in modelSBML.getListOfSpecies()]

    rxnNamesSBML = [ x.name for x in modelSBML.getListOfReactions()]
    
    for rxnID in reaction_list:
    
#         print(rxnID)
        
        rxn_info = getSpecIDs(rxnID, modelSBML, rxnNamesSBML)
        
        # rxn_params are the rows with rxnID as reaction name
        rxn_params = metabolism_params.loc[ metabolism_params["Reaction Name"] == rxnID ]

#         print(rxn_info)

    #     print(rxn_params)

        # From .xml file to read the list of substrates and products and their stoichiometry
        substrates_list = rxn_info[0][0]
        substrates_stoich = rxn_info[0][1]
        products_list = rxn_info[1][0]
        products_stoich = rxn_info[1][1]
#         print(substrates_list)
#         print(substrates_stoich)
#         print(products_list)
#         print(products_stoich)
        
        substrate_count = int(-np.sum(substrates_stoich))
        product_count = int(np.sum(products_stoich))

        # Only need the number of reactions and products to determine the ratelaw
        rateLaw = Enzymatic(substrate_count, product_count)
        
        # if rxnID == 'ADK1' or 'ADPT':
        #     print(f"rateLaw {rateLaw}")
        #     print(f"products_list {products_list}", f"substrates_list {substrates_list}")
        #     print(f"rxn_info, {rxn_info}", f"substrate_count {substrate_count}", f"product count {product_count}")
        
        rateName = rxnID+'_rate'
        
        # Do we need to addrateform per reaction?
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)

        kcatF = rxn_params.loc[ rxn_params["Parameter Type"] == "Substrate Catalytic Rate Constant" ]["Value"].values[0]
        kcatR = rxn_params.loc[ rxn_params["Parameter Type"] == "Product Catalytic Rate Constant" ]["Value"].values[0]
        
        model.addParameter(rxnIndx, 'kcatF', kcatF)
        model.addParameter(rxnIndx, 'kcatR', kcatR)
        

        rxn_KMs = rxn_params.loc[ rxn_params["Parameter Type"] == "Michaelis Menten Constant" ]
        
        # Define the reactants, products, their concentrations, and Michaelis Menten Constants
        sub_rxn_indx_counter = 0
        for i in range(len(substrates_list)):

            metID = substrates_list[i]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    # for excellular species, their concentration is fixed
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
            
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(-substrates_stoich[i])

            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
            # The stoichiometry of reactants can be more than 1
            for j in range(stoichiometry):
                
                sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
                rateFormID = 'Sub' + str(sub_rxn_indx_counter)
            
                if metID.endswith('_e'):
                    
#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                    # for excellular species, their concentration is fixed and as parameter
                    model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

                else:
#                     print(rxnIndx, rateFormID, metID)
                    if metID == 'M_o2_c':
                        model.addParameter(rxnIndx, rateFormID, metID)
                    else:
                        model.addSubstrate(rxnIndx, rateFormID, metID)
                
                KM_ID = 'KmSub' + str(sub_rxn_indx_counter)
                
                model.addParameter(rxnIndx, KM_ID, met_KM)
                
#                 print(rxnIndx, KM_ID, met_KM)




#             print(metID, stoichiometry, met_KM)
        # For product side
        prod_rxn_indx_counter = 0
        for i in range(len(products_list)):
            
            metID = products_list[i]
            
            if metID not in list(model.getMetDict().keys()):
                
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
                
                    model.addMetabolite(metID, metID, spcConc)
                
#             print('Added metabolite: ', metID, spcConc)

            stoichiometry = int(products_stoich[i])
#             print(metID, stoichiometry)

            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
            for j in range(stoichiometry):
                
                prod_rxn_indx_counter = prod_rxn_indx_counter + 1
                
                rateFormID = 'Prod' + str(prod_rxn_indx_counter)
            
                if metID.endswith('_e'):
                    
#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                    model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

                else:
#                     print(rxnIndx, rateFormID, metID)
                    model.addProduct(rxnIndx, rateFormID, metID)
#                     if j==0:
#                         if stoichiometry == 1:
#                             model.addProduct(rxnIndx, rateFormID, metID)
#                         else:
#                             print(rxnID, metID, stoichiometry)
#                             model.addProduct(rxnIndx, rateFormID, metID, stoich=int(stoichiometry))
#                     else:
#                         print(rxnID, metID, j)
#                         model.addParameter(rxnIndx, rateFormID, metID)
#                     print(j, rateFormID, metID)
                
                KM_ID = 'KmProd' + str(prod_rxn_indx_counter)
                
                model.addParameter(rxnIndx, KM_ID, met_KM)

#             print(metID, stoichiometry, met_KM)
            
        EnzymeConc = getEnzymeConc(rxn_params, sim_properties)
            
#         EnzymeConc = partTomM(rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0], sim_properties)
        
        model.addParameter(rxnIndx, "Enzyme", EnzymeConc)

        model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1)

        # Record the enzyme spec so the per-hook refresh can recompute the enzyme
        # concentration from the current counts (no Excel/SBML reread). Keyed both
        # by reaction index (no-Cython _refreshEnzymeParams) and by reaction ID
        # (Cython _buildParamVector, which maps opt-space slots -> rxnID).
        if cache is not None:
            EnzymeStr = rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0]
            GPRrule = None
            if len(EnzymeStr.split('-')) > 1:
                GPRrule = rxn_params.loc[ rxn_params["Parameter Type"] == "GPR rule" ]["Value"].values[0]
            cache['enzyme_specs'].append((rxnIndx, EnzymeStr, GPRrule))
            cache['enzyme_spec_by_rxnID'][rxnID] = (EnzymeStr, GPRrule)

        # if rxnID == 'ADK1' or 'ADPT':
        #     print(model.getReaction(rxnIndx))

        

    reptModel(model)

    return model
#########################################################################################


#########################################################################################
def defineOtherRandomBindingReactions(model, sim_properties):
    # For other random binding reactions, we give the information of metID and stoichiometry in the excel file already
    params_file = sim_properties['kinetic_params_path']
    
    RXNS_params = pd.read_excel(params_file, sheet_name='Other-Random-Binding')
    
    reaction_list = []

    for row, item in RXNS_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])
#             print(item['Reaction Name'])
            
    for rxnID in reaction_list:
    
#         print(rxnID)

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]
        
        substrate_count = int(rxn_params.loc[ rxn_params["Parameter Type"] == "Substrates" ]["Value"].values[0])
        product_count = int(rxn_params.loc[ rxn_params["Parameter Type"] == "Products" ]["Value"].values[0])
        
        rateLaw = Enzymatic(substrate_count, product_count)
        
        rateName = rxnID+'_rate'
        
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)

        kcatF = rxn_params.loc[ rxn_params["Parameter Type"] == "Substrate Catalytic Rate Constant" ]["Value"].values[0]
        kcatR = rxn_params.loc[ rxn_params["Parameter Type"] == "Product Catalytic Rate Constant" ]["Value"].values[0]
        
        model.addParameter(rxnIndx, 'kcatF', kcatF)
        model.addParameter(rxnIndx, 'kcatR', kcatR)
        

        rxn_KMs = rxn_params.loc[ rxn_params["Parameter Type"] == "Michaelis Menten Constant" ]
        
        
        for i in range(1, substrate_count+1):

            metID = rxn_params.loc[ rxn_params["Parameter Type"] == "Sub" + str(i) ]["Value"].values[0]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
                
#                 print('Added metabolite: ', metID, spcConc)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(1)
            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
#             for j in range(stoichiometry):
                
#                 sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
            rateFormID = 'Sub' + str(i)
            
            if metID.endswith('_e'):

#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

            else:
#                     print(rxnIndx, rateFormID, metID)
                model.addSubstrate(rxnIndx, rateFormID, metID)

            KM_ID = 'KmSub' + str(i)

            model.addParameter(rxnIndx, KM_ID, met_KM)
                

        for i in range(1, product_count+1):

            metID = rxn_params.loc[ rxn_params["Parameter Type"] == "Prod" + str(i) ]["Value"].values[0]
            
#             if spcID.endswith("_e"):
            
            if metID not in list(model.getMetDict().keys()):
            
                if metID.endswith('_e'):
                    
                    spcConc = sim_properties['medium'][metID]
                    
                else:
                
                    spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)
                
#                 print('Added metabolite: ', metID, spcConc)
                
                    model.addMetabolite(metID, metID, spcConc)

            stoichiometry = int(1)
            
            met_KM = rxn_KMs.loc[ rxn_KMs["Related Species"] == metID ]["Value"].values[0]
            
#             for j in range(stoichiometry):
                
#                 sub_rxn_indx_counter = sub_rxn_indx_counter + 1
                
            rateFormID = 'Prod' + str(i)
            
            if metID.endswith('_e'):

#                     spcConc = sim_properties['medium'][metID]

#                     model.addParameter(rxnIndx, rateFormID, spcConc)
                model.addParameter(rxnIndx, rateFormID, sim_properties['medium'][metID])
#                     print(rxnIndx, rateFormID, metID)

            else:
#                     print(rxnIndx, rateFormID, metID)
                model.addProduct(rxnIndx, rateFormID, metID)

            KM_ID = 'KmProd' + str(i)

            model.addParameter(rxnIndx, KM_ID, met_KM)

            
        EnzymeConc = getEnzymeConc(rxn_params, sim_properties)
            
#         EnzymeConc = partTomM(rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0], sim_properties)
        
        model.addParameter(rxnIndx, "Enzyme", EnzymeConc)
        
        model.addParameter(rxnIndx, "onoff", 1, lb=0, ub=1)
        
    return model

#########################################################################################


#########################################################################################
def defineNonRandomBindingRxns(model, sim_properties):
    # For non random binding reactions, everything is given in the excel sheet
    params_file = sim_properties['kinetic_params_path']
    
    RXNS_params = pd.read_excel(params_file, sheet_name='Non-Random-Binding Reactions')
    
    reaction_list = []
    # Get the full list of reaction names
    for row, item in RXNS_params.iterrows():
        if item['Reaction Name'] not in reaction_list:
            reaction_list.append(item['Reaction Name'])
#             print(item['Reaction Name'])
            
    for rxnID in reaction_list:
    
        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]
        # Ratelaws are given mannually
        rateLaw = str(rxn_params.loc[ rxn_params["Parameter Type"] == "Kinetic Law" ]["Value"].values[0])
        
        rateName = rxnID+'_rate'
        
        model.addRateForm(rateName, odecell.modelbuilder.RateForm(rateLaw))

        rxnIndx = model.addReaction(rxnID, rateName, rxnName="Reaction " + rxnID)
        
        for index, row in rxn_params.iterrows():
            
            param = row['Parameter Type']
            
            if (param != "Reaction Formula") and (param != "Kinetic Law"):

                if param.startswith('Sub'):
                    
                    metID = row['Value']
                    
                    if metID not in list(model.getMetDict().keys()):
            
                        if metID.endswith('_e'):

                            spcConc = sim_properties['medium'][metID]

                        else:

                            spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)

#                         print('Added metabolite: ', metID, spcConc)

                            model.addMetabolite(metID, metID, spcConc)

                    if metID.endswith('_e'):

                        model.addParameter(rxnIndx, param, sim_properties['medium'][metID])

                    else:

                        model.addSubstrate(rxnIndx, param, metID)

                elif param.startswith('Prod'):
                    
                    metID = row['Value']
                    
                    if metID not in list(model.getMetDict().keys()):
            
                        if metID.endswith('_e'):

                            spcConc = sim_properties['medium'][metID]

                        else:

                            spcConc = partTomM(sim_properties['counts'][metID][-1], sim_properties)

#                         print('Added metabolite: ', metID, spcConc)

                            model.addMetabolite(metID, metID, spcConc)

                    if metID.endswith('_e'):

                        model.addParameter(rxnIndx, param, sim_properties['medium'][metID])

                    else:

                        model.addProduct(rxnIndx, param, metID)

                else:

                    model.addParameter(rxnIndx, param, row['Value'])

    return model
#########################################################################################

# Protein not as enzymes but as reactants and products; Phosphorelay: 4 proteins; P_0227 PDH
# BaseForm in the original form of protein, listed as the first one in Metabolite IDs
# This function only defines the names of metabolites; the reactions and rates are included in others.
#########################################################################################
def addProteinMetabolites(model, sim_properties, cache=None):
    """

    Description: add the counts of each form of proteins to the ODE model
    
    For P_0621, P_0065, P_0227, their initial counts are given by initializeMetabolitesCounts
    For 4 proteins in phosphorelay, their initial counts are given by initializeProteinMetabolitesCounts

    The new generated protein from CME side will add to the base form ptsi ... in ODE
    """ 


    data_file = sim_properties['init_conc_path']
    
    ptnMets = pd.read_excel(data_file, sheet_name='Protein Metabolites')
    
    for index, row in ptnMets.iterrows():
        
        ptnID = row['Protein']
        
        metabolites = row['Metabolite IDs'].split(',')
        

        # updating the total protein counts from CME side per second
        ptnCount = sim_properties['counts'][ptnID][-1]
        
        formsCount = 0
        
        for metID in metabolites:
            
            if metID != metabolites[0]:
                
                formsCount = formsCount + sim_properties['counts'][metID][-1]
                
#                 print(metID, formsCount)
                
                model.addMetabolite(metID, metID, partTomM(sim_properties['counts'][metID][-1], sim_properties))
                
        baseFormCount = int(ptnCount - formsCount)

#         print(metabolites[0], baseFormCount)

        baseFormConc = partTomM(baseFormCount, sim_properties)

        model.addMetabolite(metabolites[0], metabolites[0], baseFormConc)

        # Record the base-form recipe so _refreshModel can recompute its count
        # (ptnCount - sum of other forms) from the current counts each hook.
        if cache is not None:
            cache['base_forms'][metabolites[0]] = (ptnID, [m for m in metabolites if m != metabolites[0]])

    return model
#########################################################################################

# Enzymatic: Define the form of Random Binding Rate Constant
#########################################################################################
def Enzymatic(subs, prods):
        
    def numerator(subs,prods):
        
        subterm = [ '( $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subNumer = ' * '.join(subterm)
        
        prodterm = [ '( $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodNumer = ' * '.join(prodterm)
        
        numerator = '( ' + '$kcatF * ' + subNumer + ' - ' + '$kcatR * ' + prodNumer + ' )'
        return numerator
    
    def denominator(subs,prods):
        
        subterm = [ '( 1 + $Sub' + str(i) + ' / $KmSub' + str(i) + ' )' for i in range(1,subs+1)]
        subDenom = ' * '.join(subterm)
        
        prodterm = [ '( 1 + $Prod' + str(i) + ' / $KmProd' + str(i) + ' )' for i in range(1,prods+1)]
        prodDenom = ' * '.join(prodterm)
        
        denominator = '( ' + subDenom + ' + ' + prodDenom + ' - 1 )'
        return denominator
        
    rate = '$onoff * $Enzyme * ( ' + numerator(subs,prods) + ' / ' + denominator(subs,prods) + ' )'
    
    return rate
#########################################################################################


#########################################################################################
def getEnzymeConc(rxn_params, sim_properties):
    """Build-path enzyme concentration: parse the spec from rxn_params (and emit
    the time-0 GPR log), then delegate the numeric computation to
    _enzymeConcFromSpec so the build and per-hook refresh share one code path."""

    EnzymeStr = rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Value"].values[0]

    Enzymes = EnzymeStr.split('-')

    GPRrule = None

    if len(Enzymes) > 1:
        # For 'or', all enzymes can catalyze the reaction (rate from the summed
        # concentrations); for 'and', the lowest-abundance enzyme sets the rate.
        GPRrule = rxn_params.loc[ rxn_params["Parameter Type"] == "GPR rule" ]["Value"].values[0]
        if sim_properties['time_second'][-1] == 0:
            rxnID = rxn_params.loc[ rxn_params["Parameter Type"] == "Eff Enzyme Count" ]["Reaction Name"].values[0]
            print(f"ODE, {rxnID} still using GPR rule {GPRrule.upper()} on proteins {','.join(_ for _ in Enzymes)}")

    return _enzymeConcFromSpec(EnzymeStr, GPRrule, sim_properties)


#########################################################################################
def _enzymeConcFromSpec(EnzymeStr, GPRrule, sim_properties):
    """Compute enzyme concentration (mM) from a cached spec + current counts.

    Mirrors getEnzymeConc's numeric logic exactly; called both at build time and
    on every per-hook refresh (see _refreshModel)."""

    Enzymes = EnzymeStr.split('-')

    if len(Enzymes) == 1: # will work for complexes like rnsBACD, potPassive, ATPSynthase

        if Enzymes[0] == 'default':
            return 0.001

        enzymeID = Enzymes[0]

        if enzymeID.startswith('P_'): # single protein P_XXXX
            Enzymecount = getEnzymeCount(sim_properties, enzymeID)
        else: # protein complex, ATPSynthase, rnsBACD importer
            Enzymecount = sim_properties['counts'][enzymeID][-1]

        return partTomM(Enzymecount, sim_properties)

    else:

        if GPRrule == 'or':
            Enzymecount = 0
            for ptnID in Enzymes:
                Enzymecount += getEnzymeCount(sim_properties, ptnID)
            return partTomM(Enzymecount, sim_properties)

        elif GPRrule == 'and':
            Enzymecounts = []
            for ptnID in Enzymes:
                Enzymecounts.append(getEnzymeCount(sim_properties, ptnID))
            Enzymecount = min(Enzymecounts)
            return partTomM(Enzymecount, sim_properties)

    print('Something went wrong getting enzyme count')

    return None

#########################################################################################
def getEnzymeCount(sim_properties, ptnID):
    """
    Input:
        sim_properties: 
        ptnID: P_XXXX

    Description: Get the count of enzymatic proteins/complexes
    """
    try:
        locusNum = ptnID.split('_')[1]
        ptn_count = 0

        if locusNum in sim_properties['locations_ptns']['peripheral membrane']: # peripheral MP

            # print(f"Peripheral Protein {locusNum} is metabolic enzyme")
            prefixs = ['P_', 'CP_']
            for prefix in prefixs:
                ptn_count += sim_properties['counts'][prefix+locusNum][-1]
        
        else:
            ptn_count = sim_properties['counts'][ptnID][-1]

    except: # not a single protein, but protein complexes
        ptn_count = sim_properties['counts'][ptnID][-1]

    return ptn_count


#########################################################################################
def getSpecIDs(rxnName, modelSBML, rxnNamesSBML):
    """
    Description: From .xml file to obtain the metID and stoichiometry information
    """
    returnList = []
    
    rxnObj = modelSBML.getReaction( rxnNamesSBML.index(rxnName) )
    
    # Use model SBML to get IDs, names, and stoichiometries for reactants
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfReactants() ]
    spcStoich = [ -1*float(x.getStoichiometry()) for x in rxnObj.getListOfReactants() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    specIDs_noH = []
    spcStoich_noH = []
    
    for i in range(len(specIDs)):
        
        metID = specIDs[i]
        
        if (metID != 'M_h_c') and (metID != 'M_h_e') and (metID != 'M_h2o_c')  and (metID != 'M_h2o_e'):
            
            specIDs_noH.append(metID)
            stoich = spcStoich[i]
            spcStoich_noH.append(stoich)
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
#     returnList.append( [specIDs, spcStoich] )
#     returnList.append( [spcNames, specIDs, spcStoich] )
    returnList.append( [specIDs_noH, spcStoich_noH] )
    
    # Now do the same for products
    specIDs = [ x.getSpecies() for x in rxnObj.getListOfProducts() ]
    spcStoich = [ float(x.getStoichiometry()) for x in rxnObj.getListOfProducts() ]
    spcNames = [ modelSBML.getSpecies( spcID ).name for spcID in specIDs]
    
    specIDs_noH = []
    spcStoich_noH = []
    
    for i in range(len(specIDs)):
        
        metID = specIDs[i]
        
        if (metID != 'M_h_c') and (metID != 'M_h_e') and (metID != 'M_h2o_c')  and (metID != 'M_h2o_e'):
            
            specIDs_noH.append(metID)
            stoich = spcStoich[i]
            spcStoich_noH.append(stoich)
    
    if np.any( np.isnan( spcStoich ) ):
        raise Exception('Invalid stoichiometry for reaction: {}'.format(rxnName)) 
    
#     returnList.append( [specIDs, spcStoich] )
#     returnList.append( [spcNames, specIDs, spcStoich] )
    returnList.append( [specIDs_noH, spcStoich_noH] )
    
    return returnList
#########################################################################################


