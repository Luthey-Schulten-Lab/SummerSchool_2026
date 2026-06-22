"""
Add Compelx Formation into CME
Author: Enguang Fu (enguang3@illinois.edu), Troy A Brier (troyb2@illinois.edu)
"""

########################  cme_complexation.py  #########################
# Builds complexation reactions for abitrary or pre-defined complexes
# Uses formatted file (csv) to read and build the complexes
# Assume that permease must be inserted into the membrane (and bind eachother) before any complexation can occur

########################################################################
import numpy as np
import pandas as pd
import rxns_CME
# from program.gen_info_rates import *
# from cme_rate_constants import *

# complexation_rate=complexation_rate*1000/cellVolume/avgdr

ptn_prefix = 'P_' #identifier used to define protein species, warning this should be the active form of the protein ex a membrane inserted protein should only have the membrane-inserted form not the cytoplasmic form (assuming you cannot begin to form the complex in the cytoplasm)
produced_prefix = 'Produced_'

########################################################################

def dimer(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):   
    produced_cmplx = produced_prefix + cmplx
    A = ptn_prefix+genes[0]
    B = ptn_prefix+genes[1]

    species = [cmplx, produced_cmplx]; counts = [count, 0]
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    sim.addReaction((A,B),(cmplx, produced_cmplx),diffusion_binding_rate)
    rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')


    return species, counts

def A2B2(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    """
    Assembly of hetero-tetramer
    Complexes Gyrase and TopoIV
    """
    produced_cmplx = produced_prefix + cmplx
    cmplx_A2 = cmplx + '_A2'; cmplx_A2B = cmplx + '_A2B'
    species = [cmplx, cmplx_A2, cmplx_A2B, produced_cmplx]; counts = [count, 0, 0, 0]
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    A = ptn_prefix+genes[0]
    B = ptn_prefix+genes[1]

    sim.addReaction((A,A),(cmplx_A2),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2, B), (cmplx_A2B), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2B, B), (cmplx, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    return species, counts

def A2B8(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    """
    Assembly of A2(B2)4
    Complex of KtrCD, Na/K Antiporter
    """
    produced_cmplx = produced_prefix+cmplx
    cmplx_A2 = cmplx + '_A2'; cmplx_B2 = cmplx + '_B2'; cmplx_A2B2 = cmplx + '_A2B2'
    cmplx_A2B4 = cmplx + '_A2B4'; cmplx_A2B6 = cmplx + '_A2B6'

    species = [cmplx, cmplx_A2, cmplx_B2, cmplx_A2B2, cmplx_A2B4, cmplx_A2B6, produced_cmplx]
    counts = [count]; counts.extend((len(species)-1)*[0])
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    A = ptn_prefix+genes[0]
    B = ptn_prefix+genes[1]

    sim.addReaction((A,A), (cmplx_A2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((B,B), (cmplx_B2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2, cmplx_B2), (cmplx_A2B2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2B2, cmplx_B2), (cmplx_A2B4), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2B4, cmplx_B2), (cmplx_A2B6), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2B6, cmplx_B2), (cmplx, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')


    return species, counts

def A2B2C(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    """
    Assembly of A2B2C
    Complex of RNDR, ribonucleotide reductase
    """
    produced_cmplx = produced_prefix + cmplx
    cmplx_A2 = cmplx + '_A2'; cmplx_B2 = cmplx + '_B2'; cmplx_A2B2 = cmplx + '_A2B2'
    species = [cmplx, cmplx_A2, cmplx_B2, cmplx_A2B2, produced_cmplx]
    counts = [count]; counts.extend((len(species)-1)*[0])
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    A = ptn_prefix+genes[0]
    B = ptn_prefix+genes[1]
    C = ptn_prefix+genes[2]

    sim.addReaction((A,A), (cmplx_A2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((B,B), (cmplx_B2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2, cmplx_B2), (cmplx_A2B2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    # sim.addReaction((cmplx_A2, cmplx_B2), (cmplx_A2B2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx_A2B2, C), (cmplx, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')


    return species, counts


def SMC_ScpAB(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    """
    Assembly of SMC_ScpAB
    Chromosome Structure Maintainace Complex
    """
    produced_cmplx = produced_prefix + cmplx
    SMC_smc2 = cmplx + '_S2' # homodimer of SMC proteins
    SMC_smc2A = cmplx + '_S2A' # SMC bound with kleisin
    SMC_B2 = cmplx + '_B2' # Dimer of ScpB

    species = [cmplx, SMC_smc2, SMC_smc2A, SMC_B2, produced_cmplx]
    counts = [count]; counts.extend((len(species)-1)*[0])
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    smc = ptn_prefix+genes[0]
    scpA = ptn_prefix+genes[1]
    scpB = ptn_prefix+genes[2]

    sim.addReaction((smc,smc), (SMC_smc2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((SMC_smc2,scpA), (SMC_smc2A), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((scpB,scpB), (SMC_B2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((SMC_smc2A, SMC_B2), (cmplx, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    return species, counts

# def dimer_activate(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties): 
#     """
#     For RNDR reaction
#     """
#     A = ptn_prefix+genes[0]
#     B = ptn_prefix+genes[1]
#     activator = ptn_prefix+genes[2]

#     cplx_silent = cmplx + '_silent'
#     species = [cplx_silent, cmplx]; counts = [0, count]
#     species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)

#     sim.defineSpecies(species)
#     sim.addReaction((A,B),(species[0]),diffusion_binding_rate)
#     rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

#     sim.addReaction((species[0], activator),(cmplx),diffusion_binding_rate)
#     rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    
#     return species, counts

def ECFModule(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    """
    
    TA: means T with A; Ap is A' subunit

    """
    T = ptn_prefix+genes[0]; A = ptn_prefix+genes[1]; Ap = ptn_prefix+genes[2]
    TA = 'ECF_TA'; TAp = 'ECF_TAp'; ECF = 'ECF'
    produced_cmplx = produced_prefix+cmplx
    species = [TA, TAp, ECF, produced_cmplx]; counts = [0, 0, count, 0]

    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    sim.addReaction((T,A),(TA),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((T,Ap),(TAp),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((TA, Ap),(ECF,produced_cmplx),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((TAp,A),(ECF,produced_cmplx),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    return species, counts

def ECF(sim,gene,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    ECF='ECF'
    S_componet = ptn_prefix+gene[0] #S,substrate domain
    
    species = [cmplx]; counts = [count]
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)

    sim.defineSpecies(species)

    cellVolume = sim_properties['volume_L'][-1]; NA = 6.022e23
    binding_rate = 1e5/(NA*cellVolume); unbinding_rate = 1e4/(NA*cellVolume)
    sim.addReaction((ECF,S_componet),(cmplx),binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((cmplx),(ECF,S_componet),unbinding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    
    return species, counts

def ABCtransporter(sim,genes,cmplx,count,cplx_path, diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    produced_cmplx = produced_prefix+cmplx

    species = [cmplx,produced_cmplx]; counts = [count,0]

    abc_df = pd.read_excel(cplx_path, sheet_name='ABC Transporter', dtype=str)

    transporter_df = abc_df[abc_df['Name'] ==  cmplx]

    suffixes = ['_T1T2', '_T1T2N1', '_T1T2N2', '_T1T2N1N2']; intermediates = [cmplx + suffix for suffix in suffixes]
    species.extend(intermediates); counts.extend([0]*len(intermediates))
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    if 'TMD1, TMD2' in transporter_df['Domain'].values: # piABC and thmppABC
        TMD12 = ptn_prefix + transporter_df[transporter_df['Domain'] == 'TMD1, TMD2']['Gene Product'].values[0]

        sim.addReaction((TMD12),(cmplx+'_T1T2'),10); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex') # fake fast reaction

    else:
        try:
            TMD1 = ptn_prefix + transporter_df[transporter_df['Domain'] == 'TMD1']['Gene Product'].values[0]
        except:
            TMD1 = ptn_prefix + transporter_df[transporter_df['Domain'] == 'TMD1, SBP']['Gene Product'].values[0] # thmppABC
        TMD2 = ptn_prefix + transporter_df[transporter_df['Domain'] == 'TMD2']['Gene Product'].values[0]

        sim.addReaction((TMD1,TMD2),(cmplx+'_T1T2'),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')


    if 'NBD1, NBD2' in transporter_df['Domain'].values: # rnsBACD
        NBD12 = ptn_prefix + transporter_df[transporter_df['Domain'] == 'NBD1, NBD2']['Gene Product'].values[0]
        sim.addReaction((cmplx+'_T1T2', NBD12), (cmplx+'_T1T2N1N2'), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    else:
        NBD1 = ptn_prefix + transporter_df[transporter_df['Domain'] == 'NBD1']['Gene Product'].values[0]
        NBD2 = ptn_prefix + transporter_df[transporter_df['Domain'] == 'NBD2']['Gene Product'].values[0]

        sim.addReaction((cmplx+'_T1T2',NBD1),(cmplx+'_T1T2N1'),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
        sim.addReaction((cmplx+'_T1T2',NBD2),(cmplx+'_T1T2N2'),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
        
        sim.addReaction((cmplx+'_T1T2N1', NBD2), (cmplx+'_T1T2N1N2'), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
        sim.addReaction((cmplx+'_T1T2N2', NBD1), (cmplx+'_T1T2N1N2'), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
        

    if 'SBP' in transporter_df['Domain'].values:
        SBP = ptn_prefix + transporter_df[transporter_df['Domain'] == 'SBP']['Gene Product'].values[0]
        sim.addReaction((cmplx+'_T1T2N1N2', SBP), (cmplx, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    else: # sprmABC
        sim.addReaction((cmplx+'_T1T2N1N2'), (cmplx, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    
    return species, counts

def SecYEGDF(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    """
    Description: Assembly of SecYEGDF
    """

    Y = ptn_prefix+genes[0]; E = ptn_prefix+genes[1]; G = ptn_prefix+genes[2]; DF = ptn_prefix+genes[3]
    YE = 'SecYE'; YEG = 'SecYEG'; YEGDF = 'SecYEGDF'
    produced_cmplx = produced_prefix+cmplx

    species = [YE,YEG, produced_cmplx]; counts = [0,0,0]
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    sim.addReaction((Y,E), (YE), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((YE, G), YEG, diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((YEG, DF), (YEGDF,produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    return species, counts

def RNAP(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    """
    Assembly of RNAP
    """

    alpha = ptn_prefix + '0645'
    beta = ptn_prefix + '0804'
    betap = ptn_prefix + '0803'
    sigma = ptn_prefix + '0407'

    RNAP_A2 = 'RNAP_A2'
    RNAP_A2B = 'RNAP_A2B'
    RNAP_A2Bp = 'RNAP_A2Bp'
    RNAP_core = 'RNAP_core'
    produced_cmplx = produced_prefix+RNAP_core

    RNAP = 'RNAP'

    species = [RNAP_A2, RNAP_A2B, RNAP_A2Bp, produced_cmplx]
    counts = [0]*len(species)
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    sim.addReaction((alpha, alpha), (RNAP_A2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((RNAP_A2, beta), (RNAP_A2B), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((RNAP_A2, betap), (RNAP_A2Bp), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((RNAP_A2B, betap), (RNAP_core, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((RNAP_A2Bp, beta), (RNAP_core, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((RNAP_core, sigma), (RNAP), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')


    return species, counts

def ATPase(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    produced_cmplx = produced_prefix + cmplx
    species = [produced_cmplx]; counts = [0]

    ATPase_prefix = 'ATPSynthase_'

    #protein locus tags are hardcoded
    beta   = ptn_prefix + '0790'
    alpha  = ptn_prefix + '0792'
    gamma  = ptn_prefix + '0791'
    delta  = ptn_prefix + '0793'
    epsilon= ptn_prefix + '0789'
    
    C=ptn_prefix + '0795'
    A=ptn_prefix + '0796'
    B=ptn_prefix + '0794'
    #FO
    B2=ATPase_prefix + 'B2'; AB2=ATPase_prefix +'AB2'
    intermediates = [B2, AB2]
    sim.addReaction((B,B), (B2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((B2,A),(AB2),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    species.extend(intermediates); counts.extend([0]*len(intermediates)) 

    ###c ring
    c_ring = ATPase_prefix + 'C10'
    intermediates=['0','1']
    for x in range(2,11):
        intermediates.append(ATPase_prefix +'C'+str(x))

    sim.addReaction((C,C),(intermediates[2]),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    for x in range(2,10):
        sim.addReaction((intermediates[x],C),(intermediates[x+1]),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
        for y in range(x,11-x):
            sim.addReaction((intermediates[x],intermediates[y]),(intermediates[x+y]),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    species.extend(intermediates[2:]); counts.extend(len(intermediates[2:])*[0])

    
    #F1
    ge=ATPase_prefix +'ge'
    a1b1=ATPase_prefix +'a1b1'
    a2b2=ATPase_prefix +'a2b2'
    a3b3=ATPase_prefix +'a3b3'
    a3b3ge=ATPase_prefix +'a3b3ge'
    a3b3geC10 = ATPase_prefix +'a3b3geC10'
    intermediates=[ge, a1b1, a2b2, a3b3, a3b3ge, a3b3geC10]

    species.extend(intermediates); counts.extend(len(intermediates)*[0])


    sim.addReaction((gamma, epsilon), (ge), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    sim.addReaction((alpha,beta),(a1b1),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((a1b1,a1b1),(a2b2),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((a2b2,a1b1),(a3b3),diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    sim.addReaction((ge,a3b3),(a3b3ge),diffusion_binding_rate) ; rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')# Atp12p
    sim.addReaction((a3b3ge, c_ring), (a3b3geC10), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    #Final combinations
    a3b3gdeC10 = ATPase_prefix +'a3b3gdeC10'
    dAB2 = ATPase_prefix +'dAB2'
    intermediates = [a3b3gdeC10, dAB2]
    species.extend(intermediates); counts.extend(len(intermediates)*[0])

    ATPase='ATPSynthase'
    species.extend([ATPase]); counts.extend([count])

    sim.addReaction((a3b3geC10, delta), (a3b3gdeC10), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((a3b3gdeC10, AB2), (ATPase, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    sim.addReaction((AB2, delta), (dAB2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((dAB2, a3b3geC10), (ATPase, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    return species, counts

def Degradosome(sim,genes,cmplx,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties):
    produced_cmplx = produced_prefix + cmplx
    species = [produced_cmplx]; counts = [0]

    Degradosome = 'Degradosome'
    
    Degradosome_prefix = 'Degradosome_'
    J1 = ptn_prefix + '0600'; J2 = ptn_prefix + '0257'
    Y = ptn_prefix + '0359'; M = ptn_prefix + '0437'

    # R = ptn_prefix + '0775'; B = ptn_prefix + '0410'
    # E = ptn_prefix + '0213'; P = ptn_prefix + '0220'

    J1J2 = Degradosome_prefix + 'J1J2'; YJ1J2 = Degradosome_prefix + 'YJ1J2'
    
    # YJ1J2MR = Degradosome_prefix + 'YJ1J2MR'
    # YJ1J2MRB = Degradosome_prefix + 'YJ1J2MRB'

    intermediates = [J1J2, YJ1J2]

    species.extend(intermediates); counts.extend([0]*len(intermediates))
    species, counts = rxns_CME.checkdefinedSpecies(sim, species, counts)
    sim.defineSpecies(species)

    # Rnase R 0775,and helicase cshB 0410 not considered

    sim.addReaction((J1, J2), (J1J2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((Y, J1J2 ), (YJ1J2), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    sim.addReaction((YJ1J2, M), (Degradosome, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')

    # sim.addReaction((YJ1J2M, E), (YJ1J2ME), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    # sim.addReaction((YJ1J2ME, P), (YJ1J2MRB), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    # sim.addReaction((YJ1J2MRB, E), (YJ1J2MRBE), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')
    # sim.addReaction((YJ1J2MRBE, P), (Degradosome, produced_cmplx), diffusion_binding_rate); rxns_CME.addReactionToMap(sim_properties, subsystem='ptnComplex')


    return species, counts


########################################################################
def build_complexation(sim,sim_properties, diffusion_binding_rate, diffusion_docking_rate):
    cplx_path =  sim_properties['cplx_path']

    cplx_df = pd.read_excel(cplx_path, sheet_name='Complexes', dtype=str)

    cplx_list = []; cplx_counts = []

    for index, row in cplx_df.iterrows():
        name = row['Name']; genes = row['Genes Products'].split(';')
        type = row['Assembly Pathway']; count = int(row['Init. Count'])

        if type == 'dimer':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = dimer(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type == 'A2B2':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = A2B2(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type == 'A2(B2)4':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = A2B8(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type == 'A2B2C':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = A2B2C(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type == 'SMC_ScpAB':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = SMC_ScpAB(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        # elif type=='dimer_activate':
        #     print('CHECKPOINT: Building',name,'complex as a dimer_activate')
        #     species, counts = dimer_activate(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
        #     cplx_list.extend(species); cplx_counts.extend(counts)
            
        elif type=='ECFModule':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = ECFModule(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type=='ECF':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = ECF(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type=='ABCTransporter':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = ABCtransporter(sim,genes,name,count,cplx_path,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type=='ATPSynthase':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = ATPase(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type=='SecYEGDF':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = SecYEGDF(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type=='RNAP':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = RNAP(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        elif type=='Degradosome':
            print('CHECKPOINT: Building',name,f'complex as a {type}')
            species, counts = Degradosome(sim,genes,name,count,diffusion_binding_rate,diffusion_docking_rate,sim_properties)
            cplx_list.extend(species); cplx_counts.extend(counts)

        else:
            if not type == 'Ribosome':
                print('WARNING: No Known Complexation Strategy for',name)
                sim.defineSpecies([name]) # Define but no assembly/formation
                cplx_list.extend([name]); cplx_counts.extend([count])


    predefined_cplx_df = pd.read_excel(cplx_path, sheet_name='Predefined Complexes', dtype=str)
    
    complexes = list(predefined_cplx_df['Name'])
    complexes_counts = [int(_) for _ in list(predefined_cplx_df['Init. Count'])]

    cplx_list.extend(complexes)
    cplx_counts.extend(complexes_counts)

    getCplxIntermediatesMap(sim_properties, cplx_list, cplx_counts)

    return cplx_list, cplx_counts


def getCplxIntermediatesMap(sim_properties, cplx_list, cplx_counts):
    """
    Description: 
        Get the protein complex intermediates/subcomplex's dictionary

    Assume initial counts of complex intermediates are all 0
    """

    cplx_dict = sim_properties['cplx_dict'] # Only Final holo Complex
    holo_cmplx = list(cplx_dict.keys())
    produced_holo_cmplx = [produced_prefix+_ for _ in holo_cmplx]

    cplx_path = sim_properties['cplx_path'] 
    transmemPtnList = sim_properties['locations_ptns']['trans-membrane']
    letter_locusNum = getLettertolocusNum()

    sub_cplx_dict = {}
    
    for i_cplx, cplx in enumerate(cplx_list):
        if cplx not in holo_cmplx+produced_holo_cmplx: # Exclude the completely assembled complex
            sub_cplx_dict[cplx] = {}

            sub_cplx_dict[cplx]['init_count'] = cplx_counts[i_cplx] # Always 0

            Stoi = mapStringtoStoi(cplx, cplx_path, letter_locusNum) # 

            sub_cplx_dict[cplx]['Stoi'] = Stoi
    

    for name, subdict in sub_cplx_dict.items():
        Stoi = subdict['Stoi']
        transMP_count = 0
        for locusNum, stoi in Stoi.items():
            if locusNum in transmemPtnList:
                transMP_count += stoi
        subdict['transMP_count'] = transMP_count
    
    sim_properties['sub_cplx_dict'] = sub_cplx_dict
    
    print(f"sub_cplx_dict {len(sub_cplx_dict)}")

    for cplx, subdict in sub_cplx_dict.items():
        print(cplx, subdict)

    return None

def getLettertolocusNum():

    letter_locusNum = {}

    ATPase_letter_locusNum = {'A':'0796', 'B':'0794', 'C':'0795',
                              'a': '0792', 'b':'0790', 'g':'0791',
                              'e':'0789', 'd': '0793'}
    KtrCD_letter_locusNum = {'A': '0685', 'B': '0686'}
    TopoIV_letter_locusNum = {'A': '0452', 'B': '0453'}
    Gyrase_letter_locusNum = {'A':'0006', 'B': '0007'}
    RNDR_letter_locusNum = {'A':'0771', 'B':'0773', 'C':'0772'}
    SMC_letter_locusNum = {'S':'0415', 'A':'0327', 'B':'0328'}
    Deg_letter_locusNum = {'Y':'0359', 'M':'0437', 'R':'0775', 'B':'0410', 'E':'0213'}

    letter_locusNum['ATPSynthase'] = ATPase_letter_locusNum; letter_locusNum['KtrCD'] = KtrCD_letter_locusNum
    letter_locusNum['TopoIV'] = TopoIV_letter_locusNum; letter_locusNum['Gyrase'] = Gyrase_letter_locusNum
    letter_locusNum['RNDR'] = RNDR_letter_locusNum; letter_locusNum['SMC'] = SMC_letter_locusNum
    letter_locusNum['Degradosome'] = Deg_letter_locusNum

    return letter_locusNum

def mapStringtoStoi(cplx, cplx_path, letter_locusNum):
    import re
    """
    Map the name string of subcomplex to Stoi dictionary
    """
    
    Stoi = {}

    if cplx.startswith('ECF_'): # ECFModule
        composition = cplx.split('_')[1]
        if composition == 'TA':
            Stoi = {'0641':1, '0642':1}
        elif composition == 'TAp':
            Stoi = {'0641':1, '0642':1}
        else:
            print(f'WARNING: ECFModule {cplx}')

    elif cplx.startswith('rnsBACD_') or cplx.find('ABC_') != -1: # ABC transporter
        # print(f"ABC Transporter {cplx}")
        abc_df = pd.read_excel(cplx_path, sheet_name='ABC Transporter', dtype=str)

        abc_transporter, composition = cplx.split('_')

        transporter_df = abc_df[abc_df['Name'] ==  abc_transporter]

        locusNums = []

        if 'TMD1, TMD2' in transporter_df['Domain'].values: # All complex intermediates of ABC transporters have T1T2
            locusNums.append( transporter_df[transporter_df['Domain'] == 'TMD1, TMD2']['Gene Product'].values[0])
        else:
            try:
                locusNums.append(transporter_df[transporter_df['Domain'] == 'TMD1']['Gene Product'].values[0])
            except:
                locusNums.append(transporter_df[transporter_df['Domain'] == 'TMD1, SBP']['Gene Product'].values[0])

            locusNums.append(transporter_df[transporter_df['Domain'] == 'TMD2']['Gene Product'].values[0])


        if composition != 'T1T2':
            if 'NBD1, NBD2' in transporter_df['Domain'].values:
                locusNums.append(transporter_df[transporter_df['Domain'] == 'NBD1, NBD2']['Gene Product'].values[0])
            else:
                if composition.find('N1') != -1:
                    locusNums.append(transporter_df[transporter_df['Domain'] == 'NBD1']['Gene Product'].values[0])
                if composition.find('N2') != -1:
                    locusNums.append(transporter_df[transporter_df['Domain'] == 'NBD2']['Gene Product'].values[0])

        for locusNum in set(locusNums):
            Stoi[locusNum] = locusNums.count(locusNum)


    elif cplx.startswith('RNAP_'): # RNAP_A2 A2B A2Bp core
        Stoi['0645'] = 2
        if cplx.endswith('B'):
            Stoi['0804'] = 1
        elif cplx.endswith('Bp'):
            Stoi['0803'] = 1
        elif cplx.endswith('core'):
            Stoi['0804'] = 1
            Stoi['0803'] = 1

    elif cplx == 'SecYE':
        Stoi = {'0652':1, '0839':1}
        
    elif cplx == 'SecYEG':
        Stoi = {'0652':1, '0839':1, '0774':1}

    else:
        for holo_cplx in letter_locusNum.keys():
            if cplx.find(holo_cplx) != -1:
                composition = cplx.split('_')[1]

                if holo_cplx == 'Degradosome': # Remove subunit J1 and J2 in degradosome since not fit extrac_.._string function
                    Stoi['0600'] = 1; Stoi['0257'] = 1
                    composition = re.sub(r'J1|J2', '', composition)

                matches = extract_numbers_from_string(composition)
                for letter, stoi_number in matches:
                    locusNum = letter_locusNum[holo_cplx][letter]
                    Stoi[locusNum] = stoi_number

    if len(Stoi) == 0:
        print(f"WARNING: Subcomplex {cplx} not matched")

    return Stoi

def extract_numbers_from_string(input_str):

    import re
    # Use regex to extract letters followed by digits, and assign 1 if no digit is present
    matches = re.findall(r'([a-zA-Z])(\d*)', input_str)

    # For each match, assign 1 if the digit part is empty
    return [(char, int(num) if num else 1) for char, num in matches]