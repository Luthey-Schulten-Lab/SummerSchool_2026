"""
Author: Enguang Fu

Date: March 2024

General functions that create genetic information processing reactions for individual genes.
Only define the species and reactions
Particle numbers are given in initializeCounts function

"""

import pandas as pd
import json
import numpy as np
import replication
import GIP_rates
import cme_complexation
import initiation as IC


#########################################################################################
def addReactionToMap(sim_properties, subsystem, locusNum = None):
    """
    Description: Adding reaction and its serial number to sim_properties['rxns_map']

    """
    if sim_properties['restartNum'][-1] == 0: # Only add reaction once
    
        rxns_number = sim_properties['rxns_numbers']
        rxns_map = sim_properties['rxns_map']

        rxns_prefix = sim_properties['rxns_prefix']
        prefix = rxns_prefix[subsystem]
        
        # get reaction name
        if subsystem in ['initiation', 'ribo_biogenesis', 'tRNACharging', 'translocation', 'ptnComplex']: # initiation_1,2,3
            rxn_name = prefix + '_' + str(rxns_number[subsystem][0])
        else:
            rxn_name = prefix + '_' + locusNum # rep_0001, rep_0002, ...

        # Avoid double adding
        if rxn_name in rxns_map.keys():
            print(f"WARNING: {rxn_name} already added to sim_properties['rxns_map']")
        else:
            total_number_rxns = sum(number_rxns[0] for number_rxns in rxns_number.values()) 

            rxns_number[subsystem][0] += 1
            rxns_map[rxn_name] = total_number_rxns      


    return None
#########################################################################################


#########################################################################################
def addGeneticInformationReactions(sim, sim_properties, gbfile):
    """
    Input: sim, sim_properties

    Return: None

    Called when restart new CME simulation

    Description: Define CME GIP species and add reactions; Create the list of species their initial counts of different GIP processes
    """
        
    genome = sim_properties['genome']

    # add chromosome replication initiation and replication reactions
    # Genes are defined in addReplication reaction.
    # Numbers of species involved in initiation are already given in this function;
    # For replication, the numbers of intergenic species are defined, numbers of genes are not.

    ini_list, ini_counts = replication.addRepInit(sim, sim_properties)

    # 493 genes replicated G_XXXX, Produced_G_XXXX, JCVISYN3A_XXXX_inter
    rep_list, rep_counts = replication.addReplication(sim,sim_properties, gbfile)
    
    # Predefine Species that will be used in Genetic Information Process
    pre_list, pre_counts = predefineSpecies(sim, sim_properties)

    currenttime_second = sim_properties['time_second'][-1]

    if currenttime_second == 0: 
        Bool_InitializeList = True
        # 493 gene all have transcription reactions RNAP_G_XXXX,  R_XXXX, and Produced_R_XXXX
        trsc_list = []; trsc_counts = []
        # 455 gene coding for proteins: P_XXXX, RNC_XXXX, Produced_P_XXXX
        translation_list = []; translation_counts = []
        # 455 gene coding for proteins: Degradosome_mRNA_XXXX, Degradated_mRNA_XXXX
        Deg_list = []; Deg_counts = []

        # to count occupied ribosome and degradosome
        sim_properties['occupied_complexes'] = {}; 
        sim_properties['occupied_complexes']['occupied_ribosome'] = int(0)
        sim_properties['occupied_complexes']['occupied_degradosome'] = int(0)
        sim_properties['occupied_complexes']['occupied_RNAP'] = int(0)
    else:
        Bool_InitializeList = False


# 'Type' = ['protein', 'ncRNA', 'gene', 'rRNA', 'tRNA', 'tmRNA']
    # Add transcription, translation, translocation and degradataion reactions
    for locusTag, locusDict in genome.items():
        locusNum = locusTag.split('_')[1]
        rnasequence = locusDict["RNAsequence"]

        if locusDict['Type'] == 'gene':
            # Nothing to three pseudo genes
            None

        elif locusDict['Type'] == 'rRNA':
            # Only Transcription for rRNAs
            species, counts = transcription(sim, sim_properties, locusNum)
            if Bool_InitializeList:
                trsc_list.extend(species); trsc_counts.extend(counts)


        elif locusDict["Type"] == 'protein':
            # Add ribosome binding and translation of each mRNA
            species, counts = transcription(sim, sim_properties, locusNum)
            if Bool_InitializeList:
                trsc_list.extend(species); trsc_counts.extend(counts)

            # Add degradosome binding and degradatation of each mRNA
            species, counts = degradation_mrna(sim, sim_properties, locusNum, rnasequence)
            if Bool_InitializeList:
                Deg_list.extend(species); Deg_counts.extend(counts)

            full_aasequence = locusDict["AAsequence"]

            species, counts = translationWithCost(sim, sim_properties, locusNum, full_aasequence)
            if Bool_InitializeList:
                translation_list.extend(species); translation_counts.extend(counts)


        elif locusDict["Type"] == 'tRNA':
            # 29 tRNAs
            # R_XXXX means free tRNA; acoording to cell 2022, we assume 40 free tRNA, 160 charged tRNA at initial conditions
            species, counts = transcription(sim, sim_properties, locusNum)
            if Bool_InitializeList:
                trsc_list.extend(species); trsc_counts.extend(counts)
            
        else:
            # tmRNA G_0158 and ncRNA G_0049 and G_0356
            species, counts = transcription(sim, sim_properties, locusNum)
            if Bool_InitializeList:
                trsc_list.extend(species); trsc_counts.extend(counts)

    # add tRNA charging, ribosome biogenesis, and protein complex formation
    tRNA_list, tRNA_counts = tRNAcharging(sim,sim_properties)
    ribo_list, ribo_counts = addRibosomeBiogenesis(sim, sim_properties, 19)
    cplx_list, cplx_counts = addPtnComplexes(sim, sim_properties, 1e5, 1e4)


    if Bool_InitializeList:
        sim_properties['pre_list'] = pre_list; sim_properties['pre_counts'] = pre_counts
        sim_properties['ini_list'] = ini_list; sim_properties['ini_counts'] = ini_counts
        sim_properties['rep_list'] = rep_list; sim_properties['rep_counts'] = rep_counts
        sim_properties['trsc_list'] = trsc_list; sim_properties['trsc_counts'] = trsc_counts
        sim_properties['translation_list'] = translation_list; sim_properties['translation_counts'] = translation_counts
        sim_properties['Deg_list'] = Deg_list; sim_properties['Deg_counts'] = Deg_counts
        sim_properties['tRNA_list'] = tRNA_list; sim_properties['tRNA_counts'] = tRNA_counts
        sim_properties['ribo_list'] = ribo_list; sim_properties['ribo_counts'] = ribo_counts
        sim_properties['cplx_list'] = cplx_list; sim_properties['cplx_counts'] = cplx_counts

        #print('The lists of GIP species and counts are passed to sim_properites Dictionary')
        
    return None

#########################################################################################

    

#########################################################################################
def transcription(sim, sim_properties, locusNum):
    
    # Already defined Species
    gene = 'G_' + locusNum    
    RNAP = 'RNAP'; RNAP_core = 'RNAP_core'; RNAP_sigma = 'P_0407'

    # Newly defined Species
    rnaID  = 'R_' + locusNum
    RNAP_gene = 'RNAP_G_' + locusNum   
    RNAProduced = 'Produced_R_' + locusNum   # The accumulative number of the newly synthesized mRNA from transcription
    
    species = [rnaID, RNAP_gene, RNAProduced]
    RNA_Init_count = getRNAInitCount(sim_properties, locusNum)
    counts = [RNA_Init_count, 0, 0]
    
    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)


    # RNAP_binding is the function to calculate the rate of binding
    # RNAP_binding is independent with metabolism
    ########################################################################
    # Calculate RNAP binding rate that is dependent on promoter strength and cellvolume
    PromoterStrength = sim_properties['promoters'][locusNum]
    # The volume of E coli is constant
    Ecoli_V = 1e-15  
    NA = 6.022e23

    cellVolume = sim_properties['volume_L'][-1]

    #binding_rate = 4*(180/765)*(PromoterStrength)*Ecoli_V*NA/11400/60
    # RNAP_binding_rate = 4*(180/765)*(PromoterStrength)*Ecoli_V*NA/11400/60/NA/cellVolume
    RNAP_binding_rate = 10*(180/765)*(PromoterStrength)*Ecoli_V*NA/11400/60/NA/cellVolume
    # 4 is the initiations per minute per rrn (ribosomal ptns) gene
    # (180/765)*(PromoterStrength) is the promoter strength of this certion gene
    # 11400: # of RNAPs, 1996 Bremer and Dennis paper

    ###############################################################################

    sim.addReaction((gene, RNAP), (RNAP_gene, RNAP_sigma), RNAP_binding_rate)
    addReactionToMap(sim_properties, subsystem='trsc_binding', locusNum= locusNum)
    #RNAP_gene = 'RP_' + locusNum + '_C1'  
    ## new_RNA_ID = 'RP_' + locusNum + '_f' + '_C1' 
    ### Already CME ##
    # Transcription rate is dependent with the counts/concentrations of NTPs in the metabolism
    # The concentrations of NTPs are stored under sim_properties['counts]

    sim.addReaction(RNAP_gene, tuple([RNAP_core,gene,rnaID, RNAProduced]),GIP_rates.TranscriptionRate(sim_properties, locusNum))
    addReactionToMap(sim_properties, subsystem='trsc_poly', locusNum= locusNum)

    sim_properties['occupied_complexes']['occupied_RNAP'] += 0

    return species, counts
#########################################################################################


#########################################################################################
def degradation_mrna(sim, sim_properties, locusNum, rnasequence): 

    # Already defined Species
    mRNA = 'R_' + locusNum
    Degradosome = 'Degradosome'
    
    # RNaseY_Dimer = 'RNaseY_Dimer'
    # RNaseJ1J2 = 'RNaseJ1J2'

    # Newly defined Species
    Deg_mRNA = 'Degradosome_mRNA_' + locusNum
    Degradated_mRNA = 'Degradated_mRNA_' + locusNum

    # RNaseY_Dimer_mRNA = 'RNaseY_Dimer_mRNA_' + locusNum
    # Degradated_mRNA = 'Degradated_mRNA_' + locusNum

    species = [Deg_mRNA, Degradated_mRNA]
    ini_degra_mRNA_count = 0
    counts = [ini_degra_mRNA_count, 0]
    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]
    Ecoli_V = 1e-15
    NA = 6.022e23

    # deg_bind_rate = 11*NA*Ecoli_V/60/2400/NA/cellVolume #2400 #7800 #1/M/s 0.00228
    deg_bind_rate = 11*NA*Ecoli_V/60/7800/NA/cellVolume #2400 #7800 #1/M/s 0.00228
    # 11 # of cleavage events of RNase E in E.coli per minute per RNase E; 
    # 2400 and 7800 are the number of mRNA in slow and fast growing E.coli 
    # deg_bind_rate is the event that mRNA binds with degradosome per RNase E per mRNA per second in cellVolume
    # rate of binding with degradosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, Degradosome),Deg_mRNA, deg_bind_rate)
    addReactionToMap(sim_properties, subsystem='deg_binding', locusNum= locusNum)

    # sim_properties['counts']['DM_' + locusNum] = 0  for communications between methods
    # mRNA degradation rate depends only on its length
    rna_deg_rate = GIP_rates.mrnaDegradationRate(sim_properties, locusNum)
    sim.addReaction(Deg_mRNA, tuple([Degradosome, Degradated_mRNA]), rna_deg_rate)
    addReactionToMap(sim_properties, subsystem='deg_depoly', locusNum= locusNum)

    sim_properties['occupied_complexes']['occupied_degradosome'] += ini_degra_mRNA_count

    return species, counts
#########################################################################################




#########################################################################################
def translationWithCost(sim, sim_properties, locusNum, aasequence):
    """
    
    Description: translation and translocation reactions for each mRNA
    """

    if locusNum in sim_properties['locations_ptns']['trans-membrane']:
        if locusNum in sim_properties['sole_YidC_ptns']:
            species, counts = YidC_translate(sim, sim_properties, locusNum, aasequence)
        else:
            species, counts = SRP_translate(sim, sim_properties, locusNum, aasequence)

    elif locusNum in sim_properties['locations_ptns']['lipoprotein']:
        species, counts = lipo_translate(sim, sim_properties, locusNum, aasequence)

    elif locusNum in sim_properties['locations_ptns']['extracellular']:
        species, counts = secreted_translate(sim, sim_properties, locusNum, aasequence)

    elif locusNum in sim_properties['locations_ptns']['peripheral membrane']:
        if locusNum in sim_properties['peri_MP_in_cyto']:
            species, counts = cyto_translate(sim, sim_properties, locusNum, aasequence)

        else:
            species, counts = peri_translate(sim, sim_properties, locusNum, aasequence)

    elif locusNum in sim_properties['locations_ptns']['cytoplasm']:
        species, counts = cyto_translate(sim, sim_properties, locusNum, aasequence)

    elif locusNum in sim_properties['locations_ptns']['unidentified']:
        species, counts = cyto_translate(sim, sim_properties, locusNum, aasequence)
        print(f"Protein {'P_'+locusNum} with unindentified localization is assumed to be translated only in cytoplasm")

    else:
        print(f"WARNING: No translation/translocation pathway for {'P_'+locusNum} ")

    return species, counts


def YidC_translate(sim, sim_properties, locusNum, aasequence):
    """
    YidC only posttranslational translocation
    """
    
    # Already defined species
    YidC = 'P_0908'
    Ribosome = 'Ribosome'
    mRNA = 'R_'+locusNum

    # Newly defined species
    CP = 'CP_' + locusNum
    MP = 'P_' + locusNum
    RNC = 'RNC_'+locusNum
    CP_Produced = 'Produced_CP_' + locusNum
    PtnProduced = 'Produced_P_' + locusNum
    YidC_CPtn = 'YidC_CP_'+locusNum

    # Counts
    locusTag = 'JCVISYN3A_'+locusNum
    species = [CP, MP, RNC, CP_Produced, PtnProduced, YidC_CPtn]
    PtnIniCount = sim_properties['proteomics_count'][locusTag]
    corrected_PtnInitCount, total_PtnInitCount, complex = correctInitPtnCount(PtnIniCount, locusNum, sim_properties)
    ini_RNC_count = 0
    counts = [0, corrected_PtnInitCount, ini_RNC_count, 0, 0, 0] # initially all in the membrane

    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]
    Ecoli_V = 1e-15
    NA = 6.022e23

    # mRNA bind with ribosome
    ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
    # 30-60 binding events per minute # 6800 is the number of ribosomes from 1996's Bremer and Dennis paper
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs
    sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)
    addReactionToMap(sim_properties, subsystem='tran_binding', locusNum= locusNum)
    
    # Translation in Cytoplasm
    full_aasequence = aasequence
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, full_aasequence)
    translationProduct = [mRNA, Ribosome, CP, CP_Produced]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    # append LEU_cost_unpaid etc to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])

    sim.addReaction(RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='tran_poly', locusNum= locusNum)

    # Binding of Cytosolic IMP with YidC
    CP_bind = 1e5/(NA*cellVolume)
    sim.addReaction((CP, YidC), YidC_CPtn, CP_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # Insertion by YidC
    YidC_insertion = 1
    sim.addReaction(YidC_CPtn, (YidC, MP, PtnProduced), YidC_insertion)
    addReactionToMap(sim_properties, subsystem='translocation')

    sim_properties['occupied_complexes']['occupied_ribosome'] += ini_RNC_count

    return species, counts


def SRP_translate(sim, sim_properties, locusNum, aasequence):
    """
    Sec/SRP cotranslational Translocation of integral MP with the edge case where SRP could not bind RNC in time leading to protein degradation
    """

    # Already defined
    SRP = 'P_0360' # Ffh
    SR = 'P_0429' # FtsY
    SecYEGDF = 'SecYEGDF' 
    mRNA = 'R_'+locusNum
    Ribosome = 'Ribosome'
    FtsH = 'P_0039'

    # Newly defined
    MP = 'P_' + locusNum
    RNC = 'RNC_'+locusNum
    MP_Produced = 'Produced_P_' + locusNum
    SRP_RNC = 'SRP_RNC_' + locusNum
    SR_SRP_RNC = 'SR_SRP_RNC_'+locusNum
    Sec_SR_SRP_RNC = 'SecYEGDF_SR_SRP_RNC_' + locusNum
    SecYEGDF_RNC = 'SecYEGDF_RNC_'+locusNum
    RNC_handover = 'RNC_Handover_' + locusNum 
    RNC_long = 'RNC_Long_'+locusNum
    Failed_Targeting = 'Failed_Targeting_'+locusNum
    CP = 'CP_' + locusNum
    CP_Produced = 'Produced_CP_' + locusNum
    FtsH_CP = 'FtsH_CP_' + locusNum
    Deg_ptn = 'Degradated_P_' + locusNum

    species = [MP, RNC, MP_Produced, SRP_RNC, SR_SRP_RNC, Sec_SR_SRP_RNC, SecYEGDF_RNC, RNC_handover, RNC_long, Failed_Targeting, CP, CP_Produced, FtsH_CP, Deg_ptn]
    # Counts of newly defined species
    locusTag = 'JCVISYN3A_'+locusNum
    PtnIniCount = sim_properties['proteomics_count'][locusTag]
    corrected_PtnInitCount, total_PtnInitCount, complex = correctInitPtnCount(PtnIniCount, locusNum, sim_properties)
    ini_RNC_count = 0
    counts = [corrected_PtnInitCount, ini_RNC_count] + [0]*(len(species)-2)

    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]
    Ecoli_V = 1e-15
    NA = 6.022e23

    #### Cotranslational Translocation of IMPs ####
    # mRNA bind with ribosome
    ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs
    sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)
    addReactionToMap(sim_properties, subsystem='tran_binding', locusNum= locusNum)

    # SRP bind with RNC
    SRP_bind = 1e5/(NA*cellVolume)
    sim.addReaction((RNC, SRP), SRP_RNC, SRP_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # SR bind with SRP_RNC
    SR_bind = 1e6/(NA*cellVolume)
    sim.addReaction((SRP_RNC, SR), SR_SRP_RNC, SR_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # Sec bind with SR_SRP_RNC
    Sec_bind = 1e5/(NA*cellVolume)
    sim.addReaction((SR_SRP_RNC, SecYEGDF), (Sec_SR_SRP_RNC), Sec_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # 
    Cargo_transfer = 0.95 # per second
    sim.addReaction((Sec_SR_SRP_RNC), (SecYEGDF_RNC, SRP, SR, RNC_handover), Cargo_transfer)
    addReactionToMap(sim_properties, subsystem='translocation')

    # Translation and Insertion into Membrane
    full_aasequence = sim_properties['genome'][locusTag]['AAsequence']
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, full_aasequence)
    translationProduct = [mRNA,Ribosome,MP, MP_Produced, SecYEGDF]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    # append LEU_cost_unpaid etc to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])
    sim.addReaction(SecYEGDF_RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='tran_poly', locusNum= locusNum)

    #### Degradation of IMPs ####
    # The edge case happens to proteins longer than 100 AAs
    if len(full_aasequence) < 100:
        print("short protein", locusNum, len(full_aasequence))
    else:
        # Over elongation of polypeptide chain above 100 aas
        longest_aa = full_aasequence[0:100]
        translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, longest_aa)
        translationProduct = [RNC_long, Failed_Targeting]
        for aa in longest_aa:
            if aa != '*':
                translationProduct.append(unpaidaaCostMap[aa])
        sim.addReaction(RNC, tuple(translationProduct), translation_rate)
        addReactionToMap(sim_properties, subsystem='translocation')

        # Finish rest polypeptide chain

        rest_aa = full_aasequence[100:]
        translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, rest_aa)
        translationProduct = [CP, CP_Produced, Ribosome, mRNA]
        for aa in rest_aa:
            if aa != '*':
                translationProduct.append(unpaidaaCostMap[aa])
        sim.addReaction(RNC_long, tuple(translationProduct), translation_rate)
        addReactionToMap(sim_properties, subsystem='translocation')

        # Degradation of not translocated integral MPs
        FtsH_bind = 1e5/(NA*cellVolume)
        sim.addReaction((CP, FtsH), (FtsH_CP), FtsH_bind)
        addReactionToMap(sim_properties, subsystem='translocation')

        FtsH_degradation = 0.22/60 # s^-1
        sim.addReaction((FtsH_CP), (FtsH, Deg_ptn), FtsH_degradation)
        addReactionToMap(sim_properties, subsystem='translocation')


    sim_properties['occupied_complexes']['occupied_ribosome'] += ini_RNC_count

    return species, counts

def lipo_translate(sim, sim_properties, locusNum, aasequence):
    """
    Description: Translation, translocation, and modification of lipoproteins
    """

    # Already defined
    SecA = 'P_0095'
    SecYEGDF = 'SecYEGDF' # SecY
    mRNA = 'R_'+locusNum
    Ribosome = 'Ribosome'
    Lgt1 = 'P_0818'; Lgt2 = 'P_0820'; Lsp = 'P_0518'
    FtsH = 'P_0039'

    # Newly defined
    RNC = 'RNC_'+locusNum
    Lipo = 'P_' + locusNum
    Produced_Lipo = 'Produced_P_' + locusNum

    SecA_RNC = 'SecA_RNC_' + locusNum
    SecYEGDF_SecA_RNC = 'SecYEGDF_SecA_RNC_' + locusNum
    PreP = 'PreP_' + locusNum
    Produced_PreP = 'Produced_PreP_' + locusNum

    CP = 'CP_' + locusNum
    RNC_long = 'RNC_Long_'+locusNum
    Failed_Targeting = 'Failed_Targeting_'+locusNum
    CP_Produced = 'Produced_CP_' + locusNum
    FtsH_CP = 'FtsH_CP_' + locusNum
    Deg_ptn = 'Degradated_P_' + locusNum

    # SecA_CP = 'SecA_CP_' + locusNum
    # SecYEGDF_SecA_CP = 'SecYEGDF_SecA_CP_' + locusNum

    Lgt1_PreP = 'Lgt1_PreP_' + locusNum
    Lgt2_PreP = 'Lgt2_PreP_' + locusNum
    ProP = 'ProP_' + locusNum
    Lsp_ProP = 'Lsp_ProP_' + locusNum 

    species = [Lipo, RNC, Produced_Lipo, SecA_RNC, SecYEGDF_SecA_RNC, PreP, Produced_PreP, CP, RNC_long, Failed_Targeting, CP_Produced, FtsH_CP, Deg_ptn, Lgt1_PreP, Lgt2_PreP, ProP, Lsp_ProP]
    
    # Counts
    locusTag = 'JCVISYN3A_'+locusNum
    PtnIniCount = sim_properties['proteomics_count'][locusTag]
    corrected_PtnInitCount, total_PtnInitCount, complex = correctInitPtnCount(PtnIniCount, locusNum, sim_properties)
    ini_RNC_count = 0
    counts = [corrected_PtnInitCount, ini_RNC_count]
    counts.extend((len(species)-len(counts))*[0])

    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]
    Ecoli_V = 1e-15
    NA = 6.022e23

    #### Cotranslational Translocation of Lipoproteins ####

    # mRNA bind with ribosome
    ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs
    sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)
    addReactionToMap(sim_properties, subsystem='tran_binding', locusNum= locusNum)

    # RNC bind with SecA
    SecA_RNC_bind = 1e5/(NA*cellVolume)
    sim.addReaction((SecA, RNC), (SecA_RNC), SecA_RNC_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # SecA_RNC bind with SecYEGDF
    SecYEGDF_SecA_RNC_bind = 1e5/(NA*cellVolume)
    sim.addReaction((SecA_RNC, SecYEGDF), (SecYEGDF_SecA_RNC), SecYEGDF_SecA_RNC_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # SecYEGDF_SecA_RNC translocate prelipoprotein into membrane
    full_aasequence = sim_properties['genome'][locusTag]['AAsequence']
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, full_aasequence)
    translationProduct = [mRNA,Ribosome,PreP, Produced_PreP, SecYEGDF, SecA]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    # append LEU_cost_unpaid etc to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])
    sim.addReaction(SecYEGDF_SecA_RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='tran_poly', locusNum= locusNum)

    # RNC finishes the translation
    full_aasequence = sim_properties['genome'][locusTag]['AAsequence']
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, full_aasequence)
    translationProduct = [mRNA,Ribosome,CP]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    # append LEU_cost_unpaid etc to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])
    sim.addReaction(RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='translocation')

    #### Degradation of Lipoproteins ####
    max_aa = 150
    # The edge case happens to proteins longer than 150 AAs
    if len(full_aasequence) < 150:
        print(f"Short lipoprotein: {locusNum} {len(full_aasequence)} long less than {max_aa} and longest chain reduced to {int(len(full_aasequence)/2)} ")
        max_aa = int(len(full_aasequence)/2)

    # Over elongation of polypeptide chain above 100 aas
    longest_aa = full_aasequence[0:max_aa]
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, longest_aa)
    translationProduct = [RNC_long, Failed_Targeting]
    for aa in longest_aa:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])
    sim.addReaction(RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='translocation')

    # Finish rest polypeptide chain
    rest_aa = full_aasequence[max_aa:]
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, rest_aa)
    translationProduct = [CP, CP_Produced, Ribosome, mRNA]
    for aa in rest_aa:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])
    sim.addReaction(RNC_long, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='translocation')

    # Degradation of not translocated integral MPs
    FtsH_bind = 1e5/(NA*cellVolume)
    sim.addReaction((CP, FtsH), (FtsH_CP), FtsH_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    FtsH_degradation = 0.22/60 # s^-1
    sim.addReaction((FtsH_CP), (FtsH, Deg_ptn), FtsH_degradation)
    addReactionToMap(sim_properties, subsystem='translocation')

    #### Modification of Lipoproteins ####
    # PreP bind with Lgt
    Lgt_PreP_bind = 1e5/(NA*cellVolume)
    sim.addReaction((PreP, Lgt1), (Lgt1_PreP), Lgt_PreP_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    sim.addReaction((PreP, Lgt2), (Lgt2_PreP), Lgt_PreP_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # Lgt_PreP to ProP
    Lgt_rate = 0.0134 # /s
    sim.addReaction((Lgt1_PreP), (Lgt1, ProP), Lgt_rate)
    addReactionToMap(sim_properties, subsystem='translocation')

    sim.addReaction((Lgt2_PreP), (Lgt2, ProP), Lgt_rate)
    addReactionToMap(sim_properties, subsystem='translocation')

    # ProP bind with Lsp
    Lsp_ProP_bind = 1e5/(NA*cellVolume)
    sim.addReaction((ProP, Lsp), (Lsp_ProP), Lsp_ProP_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # Lsp_ProP to P_
    Lsp_rate = 0.033 # /s # 2 /min
    sim.addReaction((Lsp_ProP), (Lipo, Lsp, Produced_Lipo), Lsp_rate)
    addReactionToMap(sim_properties, subsystem='translocation')


    return species, counts

def secreted_translate(sim, sim_properties, locusNum, aasequence):
    """
    Description: Translation of secreted proteins
    """

    # Already defined
    SecA = 'P_0095'
    SecYEGDF = 'SecYEGDF' # SecY
    mRNA = 'R_'+locusNum
    Ribosome = 'Ribosome'

    # Newly defined
    RNC = 'RNC_'+locusNum
    ptnID = 'P_' + locusNum
    PtnProduced = 'Produced_P_' + locusNum
    SecA_RNC = 'SecA_RNC_' + locusNum
    SecYEGDF_SecA_RNC = 'SecYEGDF_SecA_RNC_' + locusNum

    species = [ptnID, PtnProduced, RNC, SecA_RNC, SecYEGDF_SecA_RNC]
    
    # Counts
    locusTag = 'JCVISYN3A_'+locusNum
    PtnIniCount = sim_properties['proteomics_count'][locusTag]
    corrected_PtnInitCount, total_PtnInitCount, complex = correctInitPtnCount(PtnIniCount, locusNum, sim_properties)
    ini_RNC_count = 0
    counts = [corrected_PtnInitCount, 0, 0, 0, 0]

    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]
    Ecoli_V = 1e-15
    NA = 6.022e23

    # mRNA bind with ribosome
    ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs
    sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)
    addReactionToMap(sim_properties, subsystem='tran_binding', locusNum= locusNum)

    # RNC bind with SecA
    SecA_RNC_bind = 1e5/(NA*cellVolume)
    sim.addReaction((SecA, RNC), (SecA_RNC), SecA_RNC_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # SecA_RNC bind with SecYEGDF
    SecYEGDF_SecA_RNC_bind = 1e5/(NA*cellVolume)
    sim.addReaction((SecA_RNC, SecYEGDF), (SecYEGDF_SecA_RNC), SecYEGDF_SecA_RNC_bind)
    addReactionToMap(sim_properties, subsystem='translocation')

    # SecYEGDF_SecA_RNC translocate secretory protein to outer membrane
    full_aasequence = sim_properties['genome'][locusTag]['AAsequence']
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, full_aasequence)
    translationProduct = [mRNA, Ribosome, ptnID, PtnProduced, SecYEGDF, SecA]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    # append LEU_cost_unpaid etc to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])
    sim.addReaction(SecYEGDF_SecA_RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='tran_poly', locusNum= locusNum)

    sim_properties['occupied_complexes']['occupied_ribosome'] += ini_RNC_count
   
    return species, counts


def cyto_translate(sim, sim_properties, locusNum, aasequence):
    """
    Description: translation reactions without translocation; unpaid aa cost defined on the product side
    """

    # Already defined
    mRNA = 'R_'+locusNum
    Ribosome = 'Ribosome'

    # Newly defined
    ptnID = 'P_' + locusNum
    RNC = 'RNC_'+locusNum
    PtnProduced = 'Produced_P_' + locusNum
    species = [ptnID, RNC, PtnProduced]

    locusTag = 'JCVISYN3A_'+locusNum
    PtnIniCount = sim_properties['proteomics_count'][locusTag]
    ini_RNC_count = 0
    corrected_PtnInitCount, total_PtnInitCount, complex = correctInitPtnCount(PtnIniCount, locusNum, sim_properties)
    counts = [corrected_PtnInitCount, ini_RNC_count, 0]

    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]
    Ecoli_V = 1e-15 # L
    NA = 6.022e23

    # mRNA bind with ribosome
    ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)
    addReactionToMap(sim_properties, subsystem='tran_binding', locusNum= locusNum)

   # translation
    full_aasequence = sim_properties['genome'][locusTag]['AAsequence']
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, full_aasequence)
    translationProduct = [mRNA,Ribosome,ptnID, PtnProduced]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    # append LEU_cost_unpaid etc to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])

    sim.addReaction(RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='tran_poly', locusNum= locusNum)
    
    sim_properties['occupied_complexes']['occupied_ribosome'] += ini_RNC_count

    return species, counts
#########################################################################################

def peri_translate(sim, sim_properties, locusNum, aasequence):
    """
    Description: translation reactions without translocation; unpaid aa cost defined on the product side
    """

    # Already defined
    mRNA = 'R_'+locusNum
    Ribosome = 'Ribosome'

    # Newly defined
    CP = 'CP_' + locusNum
    ptnID = 'P_' + locusNum
    RNC = 'RNC_'+locusNum
    PtnProduced = 'Produced_P_' + locusNum
    species = [CP, ptnID, RNC, PtnProduced]

    locusTag = 'JCVISYN3A_'+locusNum
    PtnIniCount = sim_properties['proteomics_count'][locusTag]
    ini_RNC_count = 0
    corrected_PtnInitCount, total_PtnInitCount, complex = correctInitPtnCount(PtnIniCount, locusNum, sim_properties)
    counts = [0, corrected_PtnInitCount, ini_RNC_count, 0] # initial all in peripheral region

    species, counts = checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]
    Ecoli_V = 1e-15 # L
    NA = 6.022e23

    # mRNA bind with ribosome
    ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
    # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
    # rate of binding between ribosome and mRNA is constant for all mRNAs

    sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)
    addReactionToMap(sim_properties, subsystem='tran_binding', locusNum= locusNum)

   # translation
    full_aasequence = sim_properties['genome'][locusTag]['AAsequence']
    translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum, full_aasequence)
    translationProduct = [mRNA, Ribosome, CP, PtnProduced]
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    # append LEU_cost_unpaid etc to the product side
    for aa in aasequence:
        if aa != '*':
            translationProduct.append(unpaidaaCostMap[aa])

    sim.addReaction(RNC, tuple(translationProduct), translation_rate)
    addReactionToMap(sim_properties, subsystem='tran_poly', locusNum= locusNum)
    
    # CP to P
    CP2P_rate = 10
    sim.addReaction(CP, ptnID, CP2P_rate)
    addReactionToMap(sim_properties, subsystem='translocation')

    # P to CP
    P2CP_rate = 1
    sim.addReaction(ptnID, CP, P2CP_rate)
    addReactionToMap(sim_properties, subsystem='translocation')

    sim_properties['occupied_complexes']['occupied_ribosome'] += ini_RNC_count

    return species, counts
#########################################################################################



# #########################################################################################
# def membranePtnTranslationWithCost(sim, sim_properties, locusNum, aasequence):
#     """
    
#     Description: translation reactions with translocation; unpaid aa cost defined on the product side
#     """


#     ptnID = 'P_' + locusNum
#     cyto_ptnID = 'C_P_' + locusNum
#     translocated_ptnID = 'Located_P_' + locusNum
#     mRNA = 'R_'+locusNum
#     Ribosome = 'Ribosome'
#     RNC = 'Ribosome_mRNA_'+locusNum
#     PtnProduced = 'Produced_P_' + locusNum

#     # Avoid define P_0652 twice
#     if locusNum != '0652':
#         species = [ptnID, cyto_ptnID, translocated_ptnID, RNC, PtnProduced]
#         sim.defineSpecies(species)

#     cellVolume = sim_properties['volume_L'][-1]

#     Ecoli_V = 1e-15
#     NA = 6.022e23

#     ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
#     # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
#     # rate of binding between ribosome and mRNA is constant for all mRNAs

#     sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)
#     addReactionToMap(sim_properties, subsystem='tran_binding', locusNum= locusNum)

#     translation_rate = GIP_rates.TranslationRate(sim_properties, locusNum)

#     translationProduct = [mRNA,Ribosome,ptnID, PtnProduced]
#     unpaidaaCostMap = sim_properties['unpaidaaCostMap']

#     # append LEU_cost_unpaid to the product side
#     for aa in aasequence:
#         if aa != '*':
#             translationProduct.append(unpaidaaCostMap[aa])


#     sim.addReaction(RNC, (mRNA, Ribosome, cyto_ptnID, PtnProduced), translation_rate)
#     addReactionToMap(sim_properties, subsystem='tran_poly', locusNum= locusNum)

#     translocation_rate = GIP_rates.TranslocationRate(sim_properties, locusNum)

#     # Secy insert membrane proteins into membranes   
#     secy = 'P_0652'

#     sim.addReaction((cyto_ptnID, secy), (ptnID, secy, translocated_ptnID), translocation_rate)
    
#     addReactionToMap(sim_properties, subsystem='translocation', locusNum= locusNum)

#     return None

# #########################################################################################


# #########################################################################################
# def translation(sim, sim_properties, locusNum, aasequence):
#     """
    
#     Description: translation reactions without translocation; no unpaid aa cost defined on the product side
#     """

# # Use a two-step explicit binding model 
    
#     ptnID = 'P_' + locusNum
#     mRNA = 'R_'+locusNum
#     Ribosome = 'Ribosome'
#     RNC = 'Ribosome_mRNA_'+locusNum
#     PtnProduced = 'Produced_P_' + locusNum
#     ribo_released = 'ribosome_released'
   
#     species = [ptnID, RNC, PtnProduced]
#     sim.defineSpecies(species)

#     cellVolume = sim_properties['volume_L'][-1]

#     Ecoli_V = 1e-15
#     NA = 6.022e23

#     ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
#     # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
#     # rate of binding between ribosome and mRNA is constant for all mRNAs

#     sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)

   
#     translation_rate = GIP_rates.TranslationRate(sim_properties, aasequence)
    
#     sim.addReaction(RNC, (mRNA,Ribosome,ptnID, PtnProduced,ribo_released), translation_rate)
    


#     return None




# #########################################################################################

# def membranePtnTranslation(sim, sim_properties, locusNum, aasequence):


#     ptnID = 'P_' + locusNum
#     cyto_ptnID = 'C_P_' + locusNum
#     mRNA = 'R_'+locusNum
#     Ribosome = 'Ribosome'
#     RNC = 'Ribosome_mRNA_'+locusNum
#     PtnProduced = 'Produced_P_' + locusNum
    
#      # Avoid define P_0652 twice
#     if locusNum != '0652':
#         species = [ptnID, cyto_ptnID, RNC, PtnProduced]
#         sim.defineSpecies(species)

#     cellVolume = sim_properties['volume_L'][-1]

#     Ecoli_V = 1e-15
#     NA = 6.022e23

#     ribo_bind = 60*Ecoli_V*NA/60/6800/NA/cellVolume  # 0.002925642336248077 1/s
#     # ribo_bind is the event that mRNA binds with ribosome per ribosome per mRNA per second in cellVolume
#     # rate of binding between ribosome and mRNA is constant for all mRNAs

#     sim.addReaction((mRNA, Ribosome), RNC, ribo_bind)

#     translation_rate = GIP_rates.TranslationRate(sim_properties, aasequence)
#     sim.addReaction(RNC, (mRNA, Ribosome, cyto_ptnID, PtnProduced), translation_rate)

#     translocation_rate = GIP_rates.TranslocationRate(aasequence)

#     # Secy insert membrane proteins into membranes   
#     secy = 'P_0652'

#     sim.addReaction((cyto_ptnID, secy), (ptnID, secy), translocation_rate)

#     return None

# #########################################################################################


#########################################################################################
def tRNAcharging(sim, sim_properties):
    """
    Input: sim, sim_properties

    Return: tRNA_list:list of species in tRNA_charging, tRNA themselves not included, tRNA_counts: initial numbers of species

    Called when restarting CME simulation; Put behind the initiationMetabolites

    Description: define the species (except free tRNA and proteins) in tRNA charging reactions in CME
    add reactions for all tRNA charging reactions;
    create the lists of defined species and their initial counts lists
    
    """ 
    n_rnxs_tRNACharging = 0

    countsDic = sim_properties['counts']
    
    tRNA_list = ['M_atp_c', 'M_amp_c', 'M_ppi_c']; tRNA_counts = [countsDic['M_atp_c'][-1], countsDic['M_amp_c'][-1],countsDic['M_ppi_c'][-1]]
    sim.defineSpecies(tuple(tRNA_list))
    
    # parameters for charging of 20 amino acids
    RXNS_params = pd.read_excel( sim_properties['kinetic_params_path'], sheet_name='tRNA Charging')

    # Mapping LEU to tRNAs
    # trna_map {'LEU': ['R_0070', 'R_0423', 'R_0506'],...}
    i_isoform = 0; i_aa = 0

    for tRNA_aa, rnaIDlist in sim_properties['trna_map'].items():
        
        i_aa += 1

        tRNA_XXXX = []; tRNA_counts_XXXX = []
            
        # rxnID: LEUTRS
        rxnID = tRNA_aa + 'TRS'
        aaCost_unpaid = tRNA_aa + '_cost_unpaid'
        aaCost = tRNA_aa + '_cost'

        rxn_params = RXNS_params.loc[ RXNS_params["Reaction Name"] == rxnID ]

        # aaID: M_ala_L_c
        aaID = rxn_params.loc[ rxn_params["Parameter Type"] == 'amino acid' ]["Value"].values[0]

        aacount = countsDic[aaID][-1]
    
        # P_0163 .... already defined in translations
        synthetaseID = rxn_params.loc[ rxn_params["Parameter Type"] == 'synthetase' ]["Value"].values[0]
        # P_0163_atp
        synthetaseAtpID = synthetaseID + '_atp'
        # P_p163_atp_aa
        synthetaseAaID = synthetaseAtpID + '_aa'

        tRNA_XXXX.extend([aaID, synthetaseAtpID, synthetaseAaID])
        tRNA_counts_XXXX.extend([aacount,0,0])

        sim.defineSpecies(tuple([aaID, synthetaseAtpID, synthetaseAaID]))


        sim.addReaction(tuple([synthetaseID, 'M_atp_c']), synthetaseAtpID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_atp' ]["Value"].values[0])
        n_rnxs_tRNACharging +=1
    
        sim.addReaction(tuple([synthetaseAtpID, aaID]), synthetaseAaID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_aa' ]["Value"].values[0])
        n_rnxs_tRNACharging += 1

        for rnaID in rnaIDlist:
            
            i_isoform += 1

            synthetaseTrnaID = synthetaseAaID + '_' + rnaID

            sim.addReaction(tuple([synthetaseAaID, rnaID]), synthetaseTrnaID, rxn_params.loc[ rxn_params["Parameter Type"] == 'k_tRNA' ]["Value"].values[0])
            n_rnxs_tRNACharging += 1
            chargedTrnaID = rnaID + '_ch'

            # Marker to know the number of produced charged tRNA
            Produced_chargedTrna = 'Produced_' + chargedTrnaID
            
            sim.defineSpecies(tuple([synthetaseTrnaID, chargedTrnaID, Produced_chargedTrna]))

            tRNA_XXXX.extend([synthetaseTrnaID, chargedTrnaID, Produced_chargedTrna])
            # The initial number of 29 charged tRNA is 160

            tRNA_counts_XXXX.extend([0, 160, 0])

            sim.addReaction(synthetaseTrnaID, tuple(['M_amp_c', 'M_ppi_c', synthetaseID, chargedTrnaID, Produced_chargedTrna]), rxn_params.loc[ rxn_params["Parameter Type"] == 'k_cat' ]["Value"].values[0])
            n_rnxs_tRNACharging += 1

            sim.addReaction(tuple([aaCost_unpaid, chargedTrnaID]), tuple([rnaID, aaCost]), 1e5)
            n_rnxs_tRNACharging += 1
        tRNA_list.extend(tRNA_XXXX); tRNA_counts.extend(tRNA_counts_XXXX)

    time_second = sim_properties['time_second']

    if time_second[-1] == 0:
        print('{0} tRNA isoforms charging reactions of {1} amino acids are added into the simulation'.format(i_isoform, i_aa))
    
    for i in range(n_rnxs_tRNACharging):
        addReactionToMap(sim_properties, subsystem='tRNACharging')

    return tRNA_list, tRNA_counts
#########################################################################################


#########################################################################################
def addRibosomeBiogenesis(sim, sim_properties, NSpecies):
    """
    Input: 
    NSpecies: 19 (Zane & Talia, add linear pathway), 145 (Tyler's)

    Description: Revised based on Zane's 4DWCM
    """
    n_rxns_ribo = 0

    ribotoLocus = sim_properties['riboToLocus']
    
    cellVolume = sim_properties['volume_L'][-1]

    ribo_list = []

    NA = 6.022e23 # Avogadro's
    # 30S Small Subunit
    assemblyData = json.load(open('../input_data/oneParamMulder-local_min.json'))
    
    maxImts = NSpecies # 19 #145 # cfg.maxSpecies - (len(cfg.species) - 2) # Don't count 30S or 16S
    mses = assemblyData['netmin_rmse']
    spsRemoved = assemblyData['netmin_species']
    spNames = set([ sp['name'] for sp in assemblyData['species'] if sp['name'][0] == 'R' ])
    
    # Get NSpecies intermediates with largest fluxes
    for err,sps in zip(mses,spsRemoved):
        if len(spNames) <= maxImts:
            break
        spNames.difference_update(set(sps))
    
    #Remove 16S 'R' is 16S rRNA
    intermediates_SSU = tuple(spNames - set('R'))

    def count_s(string):
        return string.count('s')
    intermediates_SSU = sorted(intermediates_SSU, key=count_s)

    sim.defineSpecies(intermediates_SSU); ribo_list.extend(intermediates_SSU)
    
    assemblyRxns = [r for r in assemblyData['reactions'] if r['intermediate'] in spNames and r['product'] in spNames ]
    
    bindingRates = { p['rate_id']: p for p in assemblyData['parameters'] if p['rate_id'] in set(r['rate_id'] for r in assemblyRxns) }
    
    produced_SSU = ['Produced_SSU']
    sim.defineSpecies(produced_SSU); ribo_list.extend(produced_SSU)
    
    for rxn in assemblyRxns:
        rptn = 'S' + rxn['protein'][1:] # Change s to S

        intermediates = [rxn['intermediate']]
        product = rxn['product']

        if product == 'Rs3s4s5s6s7s8s9s10s11s12s13s14s15s16s17s19s20':
            product = tuple(['Rs3s4s5s6s7s8s9s10s11s12s13s14s15s16s17s19s20', produced_SSU[0]])

        macro_rate = bindingRates[rxn['rate_id']]['rate']*1e6 # per molar per second

        rate = macro_rate/NA/cellVolume

        rptn_locusNum = ribotoLocus[rptn][0]

        rptn_Name = 'P_' + rptn_locusNum

        if intermediates == ['R']:
            intermediates = ['R_0069', 'R_0534']

        for intermediate in intermediates:
            sim.addReaction(tuple([intermediate, rptn_Name]), product, rate)
            n_rxns_ribo += 1
    
    # 50S Large Subunit


    intermediates_LSU = []

    assemblyRxns = pd.read_excel(sim_properties['kinetic_params_path'], sheet_name='LSU Assembly', index_col=None)
    
    for index, rxn in assemblyRxns.iterrows():
        intermediates_LSU.extend([rxn['intermediate'], rxn['product']])
    
    # Remove 23S rRNA
    intermediates_LSU = tuple(set(intermediates_LSU) - set('R'))
    def count_L(string):
        return string.count('L') + string.count('S')
    intermediates_LSU = sorted(intermediates_LSU, key = count_L)

    sim.defineSpecies(intermediates_LSU); ribo_list.extend(intermediates_LSU)
    
    produced_LSU = ['Produced_LSU']
    sim.defineSpecies(produced_LSU); ribo_list.extend(produced_LSU)

    for index, rxn in assemblyRxns.iterrows():
        
        product = rxn['product']
        if product == 'R5SL1L2L3L4L5L6L7L9L10L11L13L14L15L16L17L18L19L20L21L22L23L24L27L28L29L31L32L33L34L35L36':
            
            product = tuple(['R5SL1L2L3L4L5L6L7L9L10L11L13L14L15L16L17L18L19L20L21L22L23L24L27L28L29L31L32L33L34L35L36', produced_LSU[0]])

        substrate = rxn['substrate']
        macro_rate = rxn['rate /microM/Sec']*1e6 # per Molar per second

        # macro_rate = bindingRates.loc[ bindingRates['Protein'] == substrate ]['rate /microM/Sec '].values[0]*1e6 # per Molar per second
        
        rate = macro_rate/NA/cellVolume
        
                # in the Genbank file (consequently genome Dict and ribotoLocus), we use L7/L12
        if substrate == 'L7':
        
            substrate = 'L7/L12'
        # 5S rRNA coded by two operons
        if substrate == '5S':
            
            substrateNames = ['R_0067', 'R_0532']

        else:
            rptn_locusNums = ribotoLocus[substrate]

            substrateNames = ['P_' + locusNum for locusNum in rptn_locusNums]
        
        if rxn['intermediate'] != 'R':
        
            intermediateNames = [rxn['intermediate']]
        # 23S rRNA
        elif rxn['intermediate'] == 'R':
            
            intermediateNames = ['R_0068', 'R_0533']
        
        for intermediate in intermediateNames:
            for substrate in substrateNames:
                sim.addReaction(tuple([substrate, intermediate]), product, rate)
                n_rxns_ribo += 1

    # LSU and SSU bind to ribosome
    # Approximate rate; large enough to make sure LSU and SSU always bind
    LSSU2Ribo_rate = 1e7/NA/cellVolume
    produced_ribosome = ['Produced_Ribosome']
    sim.defineSpecies(produced_ribosome); ribo_list.extend(produced_ribosome)

    sim.addReaction(tuple(['R5SL1L2L3L4L5L6L7L9L10L11L13L14L15L16L17L18L19L20L21L22L23L24L27L28L29L31L32L33L34L35L36', 
                           'Rs3s4s5s6s7s8s9s10s11s12s13s14s15s16s17s19s20' ]),tuple(['Ribosome', produced_ribosome[0]]), LSSU2Ribo_rate)
    n_rxns_ribo += 1

    # Assume the initial counts of all intermediates (except precursor rRNAs) to be 0
    ribo_counts = [0]*len(ribo_list)

    print('CHECKPOINT: Building Ribosome complex')

    # print('Reactions in ribosome biogenesis added at time {0}'.format(sim_properties['time_second'][-1]))

    for i in range(n_rxns_ribo):
        addReactionToMap(sim_properties, subsystem='ribo_biogenesis')

    return ribo_list, ribo_counts
    
    
def addPtnComplexes(sim, sim_properties, D_B_Membrane, D_D_Cyto_Membrane):
    """
    Input:

    Description: Add Protein Complexes in metabolism or GIP

    """
    cellVolume = sim_properties['volume_L'][-1]

    NA = 6.022e23 # Avogadro's

    diffusion_binding_rate = D_B_Membrane/(NA*cellVolume)
    
    diffusion_docking_rate = D_D_Cyto_Membrane/(NA*cellVolume)

    cplx_list, cplx_counts = cme_complexation.build_complexation(sim, sim_properties, diffusion_binding_rate, diffusion_docking_rate)

    print(f"cplx_list, {len(cplx_list)}, {cplx_list}")
    print(f"cplx_counts,  {len(cplx_counts)},{cplx_counts}")

    return cplx_list, cplx_counts

def correctInitPtnCount(PtnIniCount, locusNum, sim_properties, print_flag=False):
    """
    Correct the initial protein counts for protein in complexes
    """
    
    cplx_dict = sim_properties['cplx_dict']

    complex = []

    total_PtnIniCount = 0
    
    old_PtnIniCount = PtnIniCount

    corrected_PtnIniCount = PtnIniCount

    # Assumption about the counts of free ribosomal protein pool (both SSU and LSU) is 5% of the complete assembled ribosome
    rPtn_locusNums = sim_properties['rPtn_locusNums']
    
    if locusNum in rPtn_locusNums:
        corrected_PtnIniCount = max(25, int(PtnIniCount - 500))
        complex.append('Ribosome')
        total_PtnIniCount = corrected_PtnIniCount + 500
        if print_flag:
            print(f"P_{locusNum} in Ribosome Biogenesis Initial Count Corrected from {PtnIniCount} to {corrected_PtnIniCount}")

    else:
        for cplx, subdict in cplx_dict.items():
            StoiDict = subdict['Stoi']
            cplx_init_count = subdict['init_count']

            if locusNum in StoiDict.keys():
                stoi = StoiDict[locusNum]
                corrected_PtnIniCount = max(0, PtnIniCount - stoi*cplx_init_count)
                total_PtnIniCount += stoi*cplx_init_count
                if print_flag:
                    print(f"P_{locusNum} in {cplx} Initial Count Corrected from {PtnIniCount} to {corrected_PtnIniCount}")

                complex.append(cplx)

                PtnIniCount = corrected_PtnIniCount

        total_PtnIniCount += corrected_PtnIniCount

    if print_flag:
        print(f"P_{locusNum} in {', '.join(complex)} Initial Count Corrected from {old_PtnIniCount} to {corrected_PtnIniCount}")
        
        print(f"P_{locusNum} Total Initial Count is {total_PtnIniCount}")

    return corrected_PtnIniCount, total_PtnIniCount, complex



def getRNAInitCount(sim_properties, locusNum):
    
    locusNumtoType = sim_properties['locusNumtoType']
    type = locusNumtoType[locusNum]

    if type == 'protein':
        # mRNA count for initializing mRNA counts
        mRNAcount = pd.read_excel(sim_properties['init_conc_path'], sheet_name = 'mRNA Count')

        locusTag = 'JCVISYN3A_' + locusNum
        avg_free_mRNA_count = mRNAcount[mRNAcount['LocusTag'] == locusTag]['free'].values[0]
        avg_RNC_count = mRNAcount[mRNAcount['LocusTag'] == locusTag]['ribosome'].values[0]
        # avg_degra_mRNA_count = mRNAcount[mRNAcount['LocusTag'] == locusTag]['degradosome'].values[0]
        avg_total_mRNA_count = avg_free_mRNA_count + avg_RNC_count
        
        if avg_total_mRNA_count == 0.0:
            avg_total_mRNA_count = 0.0001*2

        RNA_init_count = np.random.poisson(avg_total_mRNA_count)

    elif type == 'tRNA':
        RNA_init_count = 40
    else: #'ncRNA', 'gene', 'rRNA', 'tmRNA'
        RNA_init_count = 1

    return RNA_init_count


def predefineSpecies(sim, sim_properties):
    """
    PreDefine species that will be often used in gene expression
    """

    # # define secy protein that is used translocation reactions for membrane proteins
    # secSRP = [ 'Ribosome_mRNA_', 'SRP_Ribosome_mRNA_', 'SR_SRP_Ribosome_mRNA_','Sec_Ribosome_mRNA_', 'P_','Produced_P_', 'RNC_Handover_' ]
    
    # secY = [prefix + '0652' for prefix in secSRP]
    
    # YidC = [prefix + '0908' for prefix in secSRP]

    # cyto = ['P_', 'Ribosome_mRNA_', 'Produced_P_']

    # ffh = [prefix + '0360' for prefix in cyto]

    # Ftsy = [prefix + '0429' for prefix in cyto]


    # RNase_R = 'P_0775'; YhaM = 'P_0437'
    counts = []
    predefined_species = []

    SRP = 'P_0360'; SR = 'P_0429'; YidC = 'P_0908'; RNAP_sigma = 'P_0407'; FtsH = 'P_0039'; SecA = 'P_0095'
    Lgt1 = 'P_0818'; Lgt2 = 'P_0820'; Lsp = 'P_0518'

    ptns = [SRP, SR, YidC, RNAP_sigma, FtsH, SecA, Lgt1, Lgt2, Lsp]
    
    for ptn in ptns:
        locusNum = ptn.split('_')[1]
        locusTag = 'JCVISYN3A_'+locusNum
        PtnIniCount = sim_properties['proteomics_count'][locusTag]
        corrected_PtnInitCount, total_PtnInitCount, complex = correctInitPtnCount(PtnIniCount, locusNum, sim_properties)
        predefined_species.append(ptn)
        counts.append(corrected_PtnInitCount)


    # define complexes that will be used in transcription, translation, degradation, and translocation
    cplx_path = sim_properties['cplx_path']

    predefined_cplx_df = pd.read_excel(cplx_path, sheet_name='Predefined Complexes', dtype=str)
    # cplx_df = pd.read_excel(cplx_path, sheet_name='Complexes', dtype=str)

    complexes = list(predefined_cplx_df['Name'])
    complexes_counts = [int(_) for _ in list(predefined_cplx_df['Init. Count'])]
    
    predefined_species.extend(complexes)
    counts.extend(complexes_counts)

    sim.defineSpecies(predefined_species)

    # define LEU_cost, LEU_cost_unpaid, etc in CME simulation beforehand, since translation and tRNA charging will use them as reactants
    unpaidaaCostMap = sim_properties['unpaidaaCostMap']
    aaCostMap = sim_properties['aaCostMap']

    for unpaidaaCost in unpaidaaCostMap.values():
        sim.defineSpecies([unpaidaaCost])
        predefined_species.append(unpaidaaCost)
        counts.append(0)

    for aaCost in aaCostMap.values():
        sim.defineSpecies([aaCost])
        predefined_species.append(aaCost)
        counts.append(0)

    return predefined_species, counts


def checkdefinedSpecies(sim, species, counts):
    """

    Description: Compare the species list with the already defined species in sim.particleMaps; 
                If defined, exclude from the species list and the correponding count

    """
    if len(species) != len(counts):
        raise ValueError('species and counts should have same length')

    defined_species = list(sim.particleMap.keys())

    for i_specie, specie in enumerate(species):
        if specie in defined_species:
            print(f"Specie already defined: {specie}")
            del species[i_specie]; del counts[i_specie]

    return species, counts