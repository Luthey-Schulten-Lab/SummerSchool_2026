"""
Author: Enguang Fu

Date: March 2024

updates CME and ODE states, calculate costs and update membrane
"""


import numpy as np
import h5py
import initiation as IC



def updateCMEcountsFile(sim, sim_properties, filename):
    """
    Input: sim_properties dictionary; filename of lm file

    Return None
      
    Called immediately after one CME run; Works in restart CME per minute

    Description: Called at the end of one CME simulation to Append the newest counts of CME species to CountsDic
    Read in the finished LM file using h5py; 
    
    """


    f=h5py.File(filename, "r")
    # Take the last second's counts
    data=f['Simulations'][str(1).zfill(7)]['SpeciesCounts'][()].transpose()[:,-1]
    
    # PP.closeLMFile(f)
    
    f.close()
    
    countsDic = sim_properties['counts']
    
    if len(data) == len(sim.particleMap):
        print(f"CHECKPOINT: Species Numbers in HDF5 file and in sim.particleMap are both {len(data)} ")
    else:
        print(f"Error: Species Numbers in HDF5 file and in sim.particleMap are {len(data)}, and {len(sim.particleMap)} ")


    for i in range(len(data)):
        speciesName = list(sim.particleMap.keys())[i]
        countsDic[speciesName].append(int(data[i]))

    return None

def updateCMEcountsHook(sim, CME_count_array, sim_properties):

    """
    
    Description: Using David's IET hook, pass the latest counts of CME species to dictionary
    """

    countsDic = sim_properties['counts']

    for specie_id in sim.particleMap.keys():
        if specie_id == 'M_atp_c':
            print('The count of {0} before CME is {1}'.format(specie_id, countsDic[specie_id][-1]))
            print('The count of {0} after CME before monomer cost is {1}'.format(specie_id, CME_count_array[specie_id]))
            print('The cost of ATP in tRNA Charging is {0}'.format(countsDic[specie_id][-1] - CME_count_array[specie_id]))

        countsDic[specie_id].append(int(CME_count_array[specie_id]))

    return None

def checkbeforeODE(sim_properties):
    """
    Description: Check the negative CME counts before ODE
    """

    countsDic = sim_properties['counts']

    for specie_id, count in countsDic.items():
        if count[-1] < 0 and not specie_id.endswith('_shortage'):
            print(f"WARNING: Before ODE Specie {specie_id} count is {count[-1]} at time {sim_properties['time_second'][-1]}")

    return None

def updateODEcounts(sim_properties, odeResults, odemodel):
    """
    Input: 
        sim_properties: Dictionary contains the counts trajectories
        odeResults: The concentration trajectories of metabolites
        odemodel: Constructed ode model by odecell  

    Return:
        None

    Called by:
        hookSimulation
        
    Description: Update the counts of metabolites after ODE simulation into the sim_properties dictionary
    """
    
    resFinal = odeResults[-1,:]

    metIDDict = odemodel.getMetDict()

    # {'M_ACP_c': 0, 'M_ACP_R_c': 1, 'M_apoACP_c': 2, ...}

    countsDic = sim_properties['counts']

    # Species in appendedList have one more element compared to other metaboites because they have been appened in updataCME or communicateCoststoMetabolism
    # so different way to update the counts after the ODE run.
    
    appendedList = ['M_ppi_c', 'M_pi_c', 'M_atp_c', 'M_ctp_c', 'M_utp_c', 'M_gtp_c', 'M_adp_c', 'M_cdp_c', 'M_udp_c'
        , 'M_gdp_c', 'M_datp_c', 'M_dttp_c', 'M_dctp_c', 'M_dgtp_c', 'M_amp_c', 'M_ump_c', 'M_cmp_c', 'M_gmp_c']
    aa_list = sim_properties['aa_list']
    appendedList.extend(aa_list)

    monomerlist = ['M_atp_c', 'M_utp_c', 'M_ctp_c', 'M_gtp_c','M_datp_c','M_dttp_c','M_dctp_c','M_dgtp_c']

    for metID, i_met in metIDDict.items():
        if resFinal[i_met] < 0:
            # print(f"WARNING: After ODE Metabolite {metID} conc is {resFinal[i_met]:.3e} mM at time {sim_properties['time_second'][-1]} second")
            None
        if not metID.endswith('_e'): # intracellular
            if metID in appendedList:

                if metID in monomerlist: # Second Pay Cost after ODE
                    payAfterODE(metID, i_met, resFinal, countsDic, sim_properties)

                else:
                    countsDic[metID][-1] = IC.mMtoPart(resFinal[i_met], sim_properties)
            else:
                countsDic[metID].append(IC.mMtoPart(resFinal[i_met], sim_properties))


    return None 

def payAfterODE(metID, i_met, resFinal, countsDic, sim_properties):
    """
    Pay the shortage of (d)NTPs monomers after the regeneration in ODE
    """

    monomer_shortage = metID+'_shortage'
    NA = 6.022e23; cellVolume = sim_properties['volume_L'][-1]
    particleTomM = 1000/(NA*cellVolume)

    # print(f"{metID}, {resFinal[i_met]} conc mM")

    count_afterODE = IC.mMtoPart(resFinal[i_met], sim_properties)

    monomer_shortage_count = countsDic[monomer_shortage][-1]
    
    # print(f"{metID}, {count_afterODE} count")

    # print(f"{monomer_shortage_count}")

    if count_afterODE > monomer_shortage_count:
        if monomer_shortage_count > 0:
            print(f"Shortage Paid: {metID} shortage {monomer_shortage_count*particleTomM} mM is paid at time {sim_properties['time_second'][-1]} second")
        count_aftersecondpay = int(count_afterODE - monomer_shortage_count)
        monomer_shortage_aftersecondpay = 0

    else:
        monomer_shortage_aftersecondpay = int(monomer_shortage_count - count_afterODE)
        count_aftersecondpay = 0
        print(f"WARNING: After Second Paying Cost, Monomer {metID} shortage is {monomer_shortage_aftersecondpay*particleTomM} mM at time {sim_properties['time_second'][-1]} second")

    # Update monomer count
    countsDic[metID][-1] = count_aftersecondpay # M_atp_c, ...
    # Update shortage_afterPay
    countsDic[monomer_shortage][-1] = monomer_shortage_aftersecondpay # # M_atp_c_shortage, ...

    return None



    
def updateODEtoCME(sim, CMECounts, sim_properties):
    """
    Input:

    Return: None

    Called by hookSimulation
    
    Description: After the ODE run, update certain CME species's counts (tRNA charging related species including 20 aas, atp, amp, and ppi) to the following CME simulation
    """

    countsDic = sim_properties['counts']

    # tRNA_ODE_list is the list of species that both in CME and ODE
    tRNA_ODE_list =  ['M_atp_c', 'M_amp_c', 'M_ppi_c']
    aa_list = sim_properties['aa_list']
    tRNA_ODE_list.extend(aa_list)

    for specie_id in tRNA_ODE_list:

        CMECounts[specie_id] = countsDic[specie_id][-1]

        if specie_id == 'M_atp_c':
            print('The count of {0} is set to {1} in CME after ODE run'.format(specie_id, CMECounts[specie_id]) )

    return None


def updateFluxes(sim_properties, odemodel, flux_avg, flux_end, boolean_end):
    """
    
    Description: save the fluxes into sim_properties['fluxes']
    """

    fluxesDict = sim_properties['fluxes']

    rxns = odemodel.getRxnList()

    for i_rxn, rxn in enumerate(rxns):
        rxn_id = rxn.getID()
        fluxesDict['F_' + rxn_id].append(flux_avg[i_rxn])

        if boolean_end:
            fluxesDict['F_' + rxn_id + '_end'].append(flux_end[i_rxn])

    
    return None


def calculateCosts(sim_properties, gbfile, locusNumtoIndex):
    """
    Input: None

    Return: None

    Called per communication step

    Description: calculate the costs of replication, transcription, translaion, degradation and translocation
    
    """

    calculateReplicationCosts(sim_properties, gbfile, locusNumtoIndex)

    calculateTranscriptionCosts(sim_properties)

    calculateDegradationCosts(sim_properties)

    calculateTranslationCosts(sim_properties)

    calculateTranslocationCosts(sim_properties)

    calculatetRNAchargingCosts(sim_properties)

    calculatePtnDegCosts(sim_properties)
    return None


def calculateReplicationCosts(sim_properties, gbfile, locusNumtoIndex):
    """
    Input: sim_properties dictionary; genebank file; locusNumtoIndex dictionary

    Output: None

    Called by hookSimulation

    Description: Calcualte the dNTPs and ATP costs in the 1 second CME simulation and append the costs to corresponding cost species in sim_properties  
    
    1 ATP hydrolysis per bp to unwind the dsDNA
    """

    # Calculate the dNTP consumption based on the RNAsequence
    # The intergenic region not included; should be improved
    dna = gbfile

    genome = sim_properties['genome']

    countsDic = sim_properties['counts']

    nuc_repcost = {'dATP_DNArep_cost':'M_datp_c', 'dTTP_DNArep_cost':'M_dttp_c', 'dCTP_DNArep_cost':'M_dctp_c', 'dGTP_DNArep_cost':'M_dgtp_c'}
    # pseudoGenes = ['JCVISYN3A_0051', 'JCVISYN3A_0546','JCVISYN3A_0602']
    pseudoGenes = sim_properties['pseudoGenes']
    ATP_replication_cost = 0
    dATP_rep_cost = 0
    dTTP_rep_cost = 0
    dCTP_rep_cost = 0
    dGTP_rep_cost = 0

    for locusTag, locusDic in genome.items():
        # Three pseudo Genes are not replicated
        if locusTag not in pseudoGenes:
            locusNum = locusTag.split('_')[1]

            Produced_gene = 'Produced_G_' +locusNum
            geneGenerated = int(countsDic[Produced_gene][-1] - countsDic[Produced_gene][-2])

                # read the index of the start and end point of the genic and intergenic region from prebuilt dictionary
            index = locusNumtoIndex[locusNum]
            # print(index)
            dnasequence = str(dna.seq[index[0]:index[1]])
                
            dATP_count = dnasequence.count('A')
            dTTP_count = dnasequence.count('T')
            dCTP_count = dnasequence.count('C')
            dGTP_count = dnasequence.count('G')
            
            # The dNTPs as monomers are used on both leading and lagging strand
            dATP_rep_cost += (dATP_count+dTTP_count)*geneGenerated
            dTTP_rep_cost += (dTTP_count+dATP_count)*geneGenerated
            dCTP_rep_cost += (dCTP_count+dGTP_count)*geneGenerated
            dGTP_rep_cost += (dGTP_count+dCTP_count)*geneGenerated

            ATP_replication_cost += (dATP_count + dTTP_count + dCTP_count + dGTP_count)*geneGenerated



    countsDic['dATP_DNArep_cost'].append(dATP_rep_cost)
    countsDic['dTTP_DNArep_cost'].append(dTTP_rep_cost)
    countsDic['dCTP_DNArep_cost'].append(dCTP_rep_cost)
    countsDic['dGTP_DNArep_cost'].append(dGTP_rep_cost)
    countsDic['ATP_DNArep_cost'].append(ATP_replication_cost)

    return None

def calculateTranscriptionCosts(sim_properties):
    """
    Input: sim_properties dictionary
    
    Return: None

    Called by hookSimulation

    Description: calculate the ATP and NTPs costs in the transcription for mRNA, tRNA and rRNA
    and append the costs to corresponding species

    1 ATP hydrolysis per bp to unwind the dsDNA
    
    """
    countsDic = sim_properties['counts']
    genome = sim_properties['genome']
    ATP_trsc_energy_cost = 0 
    ATP_mRNA_cost = 0; ATP_rRNA_cost = 0; ATP_tRNA_cost = 0
    UTP_mRNA_cost = 0; UTP_rRNA_cost = 0; UTP_tRNA_cost = 0
    CTP_mRNA_cost = 0; CTP_rRNA_cost = 0; CTP_tRNA_cost = 0
    GTP_mRNA_cost = 0; GTP_rRNA_cost = 0; GTP_tRNA_cost = 0

    nuc_trsccosts = {'ATP_mRNA_cost':'M_atp_c', 'CTP_mRNA_cost':'M_ctp_c', 'UTP_mRNA_cost':'M_utp_c', 'GTP_mRNA_cost':'M_gtp_c'}

    for locusTag, locusDic in genome.items():
        if locusDic['Type'] == 'protein':
            locusNum = locusTag.split('_')[1]
            Produced_RNA = 'Produced_R_' +locusNum
            RNAGenerated = int(countsDic[Produced_RNA][-1] - countsDic[Produced_RNA][-2])


            rnasequence = locusDic['RNAsequence']
            ATP_count = rnasequence.count('A')
            UTP_count = rnasequence.count('U')
            CTP_count = rnasequence.count('C')
            GTP_count = rnasequence.count('G')

            ATP_mRNA_cost += ATP_count*RNAGenerated
            UTP_mRNA_cost += UTP_count*RNAGenerated
            CTP_mRNA_cost += CTP_count*RNAGenerated
            GTP_mRNA_cost += GTP_count*RNAGenerated

            ATP_trsc_energy_cost += (ATP_count + UTP_count + GTP_count + CTP_count)*RNAGenerated

        elif locusDic['Type'] == 'rRNA':
            Produced_RNA = 'Produced_R_' +locusNum
            RNAGenerated = int(countsDic[Produced_RNA][-1] - countsDic[Produced_RNA][-2])


            rnasequence = locusDic['RNAsequence']
            ATP_count = rnasequence.count('A')
            UTP_count = rnasequence.count('U')
            CTP_count = rnasequence.count('C')
            GTP_count = rnasequence.count('G')

            ATP_rRNA_cost += ATP_count*RNAGenerated
            UTP_rRNA_cost += UTP_count*RNAGenerated
            CTP_rRNA_cost += CTP_count*RNAGenerated
            GTP_rRNA_cost += GTP_count*RNAGenerated
            ATP_trsc_energy_cost += (ATP_count + UTP_count + GTP_count + CTP_count)*RNAGenerated

        elif locusDic['Type'] == 'tRNA':
            Produced_RNA = 'Produced_R_' +locusNum
            RNAGenerated = int(countsDic[Produced_RNA][-1] - countsDic[Produced_RNA][-2])


            rnasequence = locusDic['RNAsequence']
            ATP_count = rnasequence.count('A')
            UTP_count = rnasequence.count('U')
            CTP_count = rnasequence.count('C')
            GTP_count = rnasequence.count('G')

            ATP_tRNA_cost += ATP_count*RNAGenerated
            UTP_tRNA_cost += UTP_count*RNAGenerated
            CTP_tRNA_cost += CTP_count*RNAGenerated
            GTP_tRNA_cost += GTP_count*RNAGenerated
            ATP_trsc_energy_cost += (ATP_count + UTP_count + GTP_count + CTP_count)*RNAGenerated

    countsDic['ATP_mRNA_cost'].append(ATP_mRNA_cost);countsDic['ATP_rRNA_cost'].append(ATP_rRNA_cost);countsDic['ATP_tRNA_cost'].append(ATP_tRNA_cost)
    countsDic['UTP_mRNA_cost'].append(UTP_mRNA_cost); countsDic['UTP_rRNA_cost'].append(UTP_rRNA_cost); countsDic['UTP_tRNA_cost'].append(UTP_tRNA_cost)
    countsDic['GTP_mRNA_cost'].append(GTP_mRNA_cost); countsDic['GTP_rRNA_cost'].append(GTP_rRNA_cost); countsDic['GTP_tRNA_cost'].append(GTP_tRNA_cost)
    countsDic['CTP_mRNA_cost'].append(CTP_mRNA_cost); countsDic['CTP_rRNA_cost'].append(CTP_rRNA_cost); countsDic['CTP_tRNA_cost'].append(CTP_tRNA_cost)

    countsDic['ATP_trsc_cost'].append(ATP_trsc_energy_cost)

    return None


def calculatetRNAchargingCosts(sim_properties):

    """
    
    
    Description: calculate the cost of ATP in CME tRNA charging 
            Directly use the difference between ATP count since the last value is directly after 1 second's CME and the -2 value is after last ODE
    """

    countsDic = sim_properties['counts']
    
    ATP_tRNACharging_cost = countsDic['M_atp_c'][-2] - countsDic['M_atp_c'][-1]

    countsDic['ATP_tRNAcharging_cost'].append(ATP_tRNACharging_cost)

    # tRNAmap = sim_properties['trna_map']

    # for tRNAaa, tRNAlist in tRNAmap.items():
    #     # For single amino acid, multiple tRNAs can act as carriers.

    #     aaCostID = tRNAaa + '_cost'
        
    #     tRNANum = len(tRNAlist)
        
    #     perchargedtRNACost = int(countsDic[aaCostID][-1]/tRNANum)

    #     for tRNA in tRNAlist:
    #         chargedtRNAID = tRNA + '_ch'
    #         if countsDic[chargedtRNAID][-1] >= perchargedtRNACost:

    #             countsDic[chargedtRNAID][-1] = countsDic[chargedtRNAID][-1] - perchargedtRNACost

    #         else:
    #             # print(chargedtRNAID + ' runs out at time ' + str(sim_properties['time_second'][-1]))
    #             countsDic[chargedtRNAID][-1] = 0


    return None


def calculateDegradationCosts(sim_properties):

    """
    Input: sim_properties Dictionary

    Return: None

    Called by hookSimulation

    Description: Calculate the cost of ATP and productions of NMPs in degradation of mRNAs and append costs into corresponding species
    
    1 ATP hydrolysis per bp to degradate the mRNA 
    
    """

    countsDic = sim_properties['counts']
    genome = sim_properties['genome']
    NMP_recycle_counters = {'AMP_mRNAdeg_recycled':'M_amp_c', 'UMP_mRNAdeg_recycled':'M_ump_c', 'CMP_mRNAdeg_recycled':'M_cmp_c', 'GMP_mRNAdeg_recycled':'M_gmp_c'}

    ATP_mRNAdeg_cost = 0

    AMP_mRNAdeg_recycled = 0
    UMP_mRNAdeg_recycled = 0
    CMP_mRNAdeg_recycled = 0
    GMP_mRNAdeg_recycled = 0

    for locusTag, locusDict in genome.items():
        
        if locusDict["Type"] == 'protein':
            locusNum = locusTag.split('_')[1]
            Degradated_mRNA = 'Degradated_mRNA_' + locusNum
            mRNADegradated = int(countsDic[Degradated_mRNA][-1]- countsDic[Degradated_mRNA][-2])

            rnasequence = locusDict['RNAsequence']
            AMP_count = rnasequence.count('A')
            UMP_count = rnasequence.count('U')
            CMP_count = rnasequence.count('C')
            GMP_count = rnasequence.count('G')


            AMP_mRNAdeg_recycled += AMP_count*mRNADegradated
            UMP_mRNAdeg_recycled += UMP_count*mRNADegradated
            CMP_mRNAdeg_recycled += CMP_count*mRNADegradated
            GMP_mRNAdeg_recycled += GMP_count*mRNADegradated

            ATP_mRNAdeg_cost += (AMP_count + UMP_count + CMP_count + GMP_count)*mRNADegradated

    
    countsDic['AMP_mRNAdeg_recycled'].append(AMP_mRNAdeg_recycled)
    countsDic['UMP_mRNAdeg_recycled'].append(UMP_mRNAdeg_recycled)
    countsDic['CMP_mRNAdeg_recycled'].append(CMP_mRNAdeg_recycled)
    countsDic['GMP_mRNAdeg_recycled'].append(GMP_mRNAdeg_recycled)

    
    countsDic['ATP_mRNAdeg_cost'].append(ATP_mRNAdeg_cost)



    return None



def calculateTranslationCosts(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Called by hookSimulation

    Description: Calculate the GTP and amino acids with charged tRNA costs of the tranlation reactions based on the number of Produced_P_XXXX 
    
    """
    countsDic = sim_properties['counts']

    genome = sim_properties['genome']
    GTP_translate_cost = 0

    # aaCostMap = sim_properties['aaCostMap']
    

    # aaCostCounts =   {'ALA_cost': 0,'ARG_cost': 0,'ASN_cost': 0,'ASP_cost': 0,'CYS_cost': 0,'GLU_cost': 0,'GLN_cost': 0,'GLY_cost': 0,'HIS_cost': 0,
    # 'ILE_cost': 0,'LEU_cost': 0,'LYS_cost': 0,'MET_cost': 0,'PHE_cost': 0,'PRO_cost': 0,'SER_cost': 0,'THR_cost': 0,'TRP_cost': 0,'TYR_cost': 0,
    # 'VAL_cost': 0}


    for locusTag, locusDict in genome.items():
        
        if locusDict["Type"] == 'protein':

            locusNum = locusTag.split('_')[1]
            Produced_Ptn = "Produced_P_" + locusNum
            proteinGenerated = int(countsDic[Produced_Ptn][-1] - countsDic[Produced_Ptn][-2])

            aasequence = locusDict["AAsequence"]

            GTP_translate_cost += proteinGenerated*(len(aasequence) - aasequence.count('*'))*2
            
            # for aa, aaCostStr in aaCostMap.items():
            #     # Pat attention to do accumulation of keys in dictionary
            #     aaCostCounts[aaCostStr] += proteinGenerated*aasequence.count(aa)

    # for aaCostStr, aaCost in aaCostCounts.items():
    #     # print(aaCostStr, aaCost)
    #     countsDic[aaCostStr].append(aaCost)

    countsDic['GTP_translate_cost'].append(GTP_translate_cost)

    return None

def calculateTranslocationCosts(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Called by hookSimulation

    Description: calculate and append the ATP and GTP cost for translocation of trans-membrane proteins and lipoproteins
    """
    countsDic = sim_properties['counts']

    # GTP cost for trans-membrane proteins
    transmemPtnList = sim_properties['locations_ptns']['trans-membrane']
    GTP_transloc_cost = 0

    for locusNum in transmemPtnList:
        if locusNum not in sim_properties['sole_YidC_ptns']: # F0 c subunit via YidC mechanism

            RNC_handover = 'RNC_Handover_' + locusNum

            GTP_transloc_cost +=  2*(int(countsDic[RNC_handover][-1] - countsDic[RNC_handover][-2]))

    countsDic['GTP_transloc_cost'].append(GTP_transloc_cost)

    #  print(f"{sim_properties['time_second'][-1]}, GTP GTP_transloc_cost {GTP_transloc_cost}")

    # ATP cost for lipoproteins and secreted proteins
    lipoproteinList = sim_properties['locations_ptns']['lipoprotein'] 

    ATP_transloc_cost = 0

    for locusNum in lipoproteinList:
        
        Produced_PreP = 'Produced_PreP_' + locusNum

        full_aasequence = sim_properties['genome']['JCVISYN3A_'+locusNum]['AAsequence']

        ATP_transloc_cost += int(len(full_aasequence)*(countsDic[Produced_PreP][-1] - countsDic[Produced_PreP][-2])/10) # 1 ATP per 10 AAs


    for locusNum in sim_properties['locations_ptns']['extracellular']:
        Produced_ptn = 'Produced_P_' + locusNum
        
        full_aasequence = sim_properties['genome']['JCVISYN3A_'+locusNum]['AAsequence']

        ATP_transloc_cost += int(len(full_aasequence)*(countsDic[Produced_ptn][-1] - countsDic[Produced_ptn][-2])/10) # 1 ATP per 10 AAs

    countsDic['ATP_transloc_cost'].append(ATP_transloc_cost)

    return None
    
def calculatePtnDegCosts(sim_properties):
    """
    Description: Calculate the ATP cost and AA recycled in protein degradation
    """
    countsDic = sim_properties['counts']
    transmemPtnList = sim_properties['locations_ptns']['trans-membrane']
    lipoproteinList = sim_properties['locations_ptns']['lipoprotein']
    IMPYidCList = sim_properties['sole_YidC_ptns']

    Ptn_degradable = (set(transmemPtnList) | set(lipoproteinList)) - set(IMPYidCList)

    recycledaaMap = sim_properties['recycledaaMap']

    recycledaaCount = {}

    for recycledaa in recycledaaMap.values():
        recycledaaCount[recycledaa] = 0

    ATP_Ptndeg_cost = 0
    
    for locusNum in Ptn_degradable:
        locusTag = 'JCVISYN3A_' + locusNum
        Deg_ptn = 'Degradated_P_' + locusNum
        full_aasequence =  sim_properties['genome'][locusTag]['AAsequence']
        ptndeg_count = int(countsDic[Deg_ptn][-1] - countsDic[Deg_ptn][-2])

        ATP_Ptndeg_cost += len(full_aasequence)*ptndeg_count

        for letter, recycledaa in recycledaaMap.items():
            aa_count = full_aasequence.count(letter)*ptndeg_count
            recycledaaCount[recycledaa] += aa_count

    for recycledaa, aa_count in recycledaaCount.items():
        countsDic[recycledaa].append(aa_count)
    
    countsDic['ATP_Ptndeg_cost'].append(ATP_Ptndeg_cost)

    return None

def communicateCostsToMetabolism(sim_properties):

    """
    Input: sim_properties dictionary

    Returns: None
    
    Called by: hookSimulation

    Description: Sum up the recorded consumpution of nucleotides in different processes and extract from the monomer pool 
    (such as ATP_trsc_cost, ATP_mRNA_cost into M_atp_c)
        Note: Cost of ATP in tRNA charging represented explicitly in CME
    """

    # Species M_ppi_c, M_pi_c, M_atp_c, M_ctp_c, M_utp_c, M_gtp_c, M_adp_c, M_cdp_c, M_udp_c, M_gdp_c, M_datp_c, M_dttp_c, M_dctp_c, M_dgtp_c,
    # M_amp_c, M_ump_c, M_cmp_c, M_gmp_c 's tragectories are appended by one new element.
    
    countsDic = sim_properties['counts']
    time_second = sim_properties['time_second'][-1]
    NA = 6.022e23; cellVolume = sim_properties['volume_L'][-1]
    particleTomM = 1000/(NA*cellVolume)

    ###################################
    #####  NTP Costs ####
    ###################################
    # the cost of atp in tRNA charging is already represented in the CME reactions
    NTPsCostMap = {'M_atp_c':{'syn_cost':['ATP_mRNA_cost', 'ATP_tRNA_cost', 'ATP_rRNA_cost'], 'energy_cost': ['ATP_trsc_cost', 'ATP_mRNAdeg_cost', 'ATP_DNArep_cost', 'ATP_Ptndeg_cost']},
                    'M_gtp_c':{'syn_cost': ['GTP_mRNA_cost','GTP_tRNA_cost','GTP_rRNA_cost'], 'energy_cost':['GTP_translate_cost', 'GTP_transloc_cost']},
                    'M_ctp_c':{'syn_cost': ['CTP_mRNA_cost','CTP_tRNA_cost','CTP_rRNA_cost']},
                    'M_utp_c':{'syn_cost': ['UTP_mRNA_cost','UTP_tRNA_cost','UTP_rRNA_cost']},
                    }

   
    # ppi is generated in synthesis reactions.
    ppiCost = 0

    # pi is generated in energy related reactions.
    piCost = 0

    for NTPID, subDic in NTPsCostMap.items():
        
        nucName = NTPID[-5]

        NDPID = 'M_' + nucName + 'dp_c'
        
        NDP_count = countsDic[NDPID][-1]

        NTP_shortage = NTPID + '_shortage'
        
        NTP_shortage_count = countsDic[NTP_shortage][-1]
        
        # Pure NTP cost within 1 second
        NTPCost_energy = 0

        NTPCost_syn = 0

        # Counts of NTP from last second
        NTP_count = countsDic[NTPID][-1]
        
        # Sum up all cost from different processes
        for function, processes in subDic.items():
            # print(function, processes)
            
            if function == 'energy_cost':
                for process in processes:
                    
                    NTPCost_energy += countsDic[process][-1]
                    piCost += countsDic[process][-1] # pi generated in energetic reactions

            else:
                for process in processes:

                    NTPCost_syn += countsDic[process][-1]                    
                    ppiCost += countsDic[process][-1] # ppi generated in polymerization reactions

        NTP_cost_total = NTPCost_syn + NTPCost_energy

        if NTP_count > NTP_cost_total + NTP_shortage_count: # We are able to pay the cost
            NTP_count_afterpay = int(NTP_count - NTP_cost_total - NTP_shortage_count)
            NTP_shortage_count_afterpay = 0
            NDP_count = NDP_count + NTPCost_energy
        
        else: # the shortage happens
            NTP_shortage_count_afterpay = int(NTP_cost_total + NTP_shortage_count - NTP_count + 1)
            NTP_count_afterpay = 1 # Monomer count is 1 and passed to ODE if shortage happens
            NDP_count = NDP_count + NTPCost_energy
            print(f"WARNING: Shortage of Monomer {NTPID} with concentration {NTP_shortage_count_afterpay*particleTomM} mM happens at {time_second} second")

        # First Pay cost before ODE
        if NTPID == 'M_atp_c':
            # Since M_atp_c has been appended in updataCMEcounts functions
            countsDic[NTPID][-1] = NTP_count_afterpay
            # print('The count of {0} after CME is {1}'.format(NTPID, countsDic[NTPID][-1]))
            # print('Other Cost of ATP is {0}'.format(NTPCost_energy + NTPCost_syn))
            countsDic[NTP_shortage].append(NTP_shortage_count_afterpay)

        else:
            countsDic[NTPID].append(NTP_count_afterpay)
            countsDic[NTP_shortage].append(NTP_shortage_count_afterpay)


        # NDP generated in energetic reactions
        countsDic[NDPID].append(NDP_count)

    ###################################
    #####  dNTP Costs ####
    ###################################

    # dNTPs only as Monomers so no pi generated 
    dNTPsCostMap = {'M_datp_c':'dATP_DNArep_cost', 'M_dttp_c': 'dTTP_DNArep_cost','M_dctp_c':'dCTP_DNArep_cost', 'M_dgtp_c':'dGTP_DNArep_cost'}
    
    for dNTPID, process in dNTPsCostMap.items():
        
        dNTP_shortage = dNTPID + '_shortage'

        dNTP_shortage_count = countsDic[dNTP_shortage][-1]

        dNTP_count = countsDic[dNTPID][-1]

        dNTP_cost = countsDic[process][-1]

        if dNTP_count > dNTP_shortage_count + dNTP_cost:
            dNTP_shortage_count_afterpay = 0
            dNTP_count_afterpay = int(dNTP_count - dNTP_shortage_count - dNTP_cost)
        else: # Shortage happens
            dNTP_shortage_count_afterpay = int(dNTP_shortage_count + dNTP_cost - dNTP_count + 1)
            dNTP_count_afterpay = 1
            print(f"WARNING: Shortage of Monomer {dNTPID} with concentration {dNTP_shortage_count_afterpay*particleTomM} mM happens at {time_second} second")

        countsDic[dNTPID].append(dNTP_count_afterpay)

        countsDic[dNTP_shortage].append(dNTP_shortage_count_afterpay)

        ppiCost += dNTP_cost
    
    ###################################
    #####  NMP Costs/generation ####
    ###################################

    AMPsCostMap = {'M_amp_c':'AMP_mRNAdeg_recycled', 'M_ump_c':'UMP_mRNAdeg_recycled', 'M_cmp_c': 'CMP_mRNAdeg_recycled','M_gmp_c': 'GMP_mRNAdeg_recycled'  }

    for AMPID, process in AMPsCostMap.items():

        AMPCount = countsDic[AMPID][-1]
        
        AMP_recycled = countsDic[process][-1]

        AMPCount = AMPCount + AMP_recycled
        
        if AMPID == 'M_amp_c': # amp generated in tRNA charging
            countsDic[AMPID][-1] = AMPCount
        else:
            countsDic[AMPID].append(AMPCount)

    
    # pi not in CME
    countsDic['M_pi_c'].append(countsDic['M_pi_c'][-1] + piCost)
    
    # ppi in CME, already appened
    countsDic['M_ppi_c'][-1] = countsDic['M_ppi_c'][-1] + ppiCost

    ###################################
    #####  Amino Acids Recycled ####
    ###################################
    recycledaaMap = sim_properties['recycledaaMap']
    lettertoAAMet = sim_properties['lettertoAAMet']

    for letter in lettertoAAMet.keys():
        AAMet = lettertoAAMet[letter]
        recycledaa = recycledaaMap[letter]
        countsDic[AAMet][-1] = countsDic[AAMet][-1] + countsDic[recycledaa][-1] # Recycled added before ODE run

    return None


#########################################################################################
def updateSA(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Description: Calculate and update the current cell membrane surface and volume
    """
    
    countsDic = sim_properties['counts']
    
#     lipidSizes = {
#     'M_clpn_c':0.4,
#     'M_chsterol_c':0.35, # 0.35, test for Chol. value smaller
#     'M_sm_c':0.45,
#     'M_pc_c':0.55,
#     'M_pg_c':0.6,
#     'M_galfur12dgr_c':0.6,
#     'M_fa_c':0.5, # Should this be here??
#     'M_12dgr_c':0.5, # Scale down, should cdp-dag be added???
#     'M_pa_c':0.5,
#     'M_cdpdag_c':0.5,
#     }
    
    lipidSizes = {
    'M_clpn_c':0.4,
    'M_chsterol_c':0.35, 
    'M_sm_c':0.45,
    'M_pc_c':0.55,
    'M_pg_c':0.6,
    'M_galfur12dgr_c':0.6,
    'M_12dgr_c':0.5, 
    'M_pa_c':0.5,
    'M_cdpdag_c':0.5,
    }
    # 0.4 nm^2
    
    lipidSA = 0
    
    for lipid, size in lipidSizes.items():
        
        lipidSA += countsDic[lipid][-1]*size
    
    # Lipid accounts for 0.513 weight of surface area.
    lipidSA = int(lipidSA*0.513)
    
    sim_properties['SA']['SA_lipid'].append(lipidSA)

    # MP
    MP_total_count, MP_cplx_count, MP_subcplx_count, MP_translocon_RNC, MP_free_count = IC.getMPCounts(sim_properties)

    avg_MP_area = sim_properties['avg_MP_area']
    
    MPSA = int(avg_MP_area*MP_total_count) # nm^2

    sim_properties['SA']['SA_ptn'].append(MPSA)
    
    sim_properties['SA']['SA_nm2'].append(int(lipidSA + MPSA))
    
    sim_properties['SA']['SA_m2'].append(int(lipidSA + MPSA)/1e18) # Unit m^2
    

    radius_2V = 2.52e-7  # m

    cyto_radius_nm_equivalent_sphere = np.sqrt(sim_properties['SA']['SA_m2'][-1]/(4*np.pi))   #m
    
    cyto_radius_nm = min(cyto_radius_nm_equivalent_sphere, radius_2V)


    sim_properties['volume_L'].append((4/3)*np.pi*(cyto_radius_nm)**3*1000)
    
    if cyto_radius_nm_equivalent_sphere > radius_2V:
        
        sim_properties['division_started'].append(True)
    
    # print('SA: ', sim_properties['SA'])
    # print('V: ', sim_properties['volume'])
    
    # print('cyto radius: ', sim_properties['cyto_radius'])
    # print('cyto radius nm: ', sim_properties['cyto_radius_nm'])
    
    return None


