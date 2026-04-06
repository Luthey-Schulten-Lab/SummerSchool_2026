"""
Description: Functions to analyze the ribosome biogenesis in WCM

"""
import numpy as np
import json
import pandas as pd

def getSSUConnections(SSUfile_dir, N_inters):
    """
    
    Description: 
    SSU_up_dict: key as the ID of the intermediates, value are the list of the upstream precursors 
    SSU_down_dict: key as the ID of the intermediates, value are the list of the downstream precursors
    """
    
    SSU_up_dict = {}; SSU_down_dict = {}

    intermediates_SSU = getSSUInters(SSUfile_dir, N_inters)

    intermediates_SSU[0] = 'R'

    smallAssemblyFile = SSUfile_dir
    assemblyData = json.load(open(smallAssemblyFile))

    assemblyRxns = [r for r in assemblyData['reactions'] if r['intermediate'] in intermediates_SSU and r['product'] in intermediates_SSU ]

    for inter in intermediates_SSU:
        
        inter_Rxns_up = [r for r in assemblyRxns if r['product'] == inter ]
        inter_Rxns_down = [r for r in assemblyRxns if r['intermediate'] == inter]

        if inter == 'R':
            inter = '16S'

        SSU_up_dict[inter] = [r['protein'] for r in inter_Rxns_up]

        SSU_down_dict[inter] = [r['protein'] for r in inter_Rxns_down]

    return SSU_up_dict, SSU_down_dict


def getLSUConnections(LSUfile_dir):
    """
    
    Description: 
    LSU_up_dict: key as the ID of the intermediates, value are the list of the upstream precursors 
    LSU_down_dict: key as the ID of the intermediates, value are the list of the downstream precursors
    """

    LSU_up_dict, LSU_down_dict = {}, {}

    intermediates_LSU = getLSUInters(LSUfile_dir)

    intermediates_LSU[0] = 'R'

    assemblyRxns = pd.read_excel(LSUfile_dir, sheet_name='LSU reactions', index_col=None)

    for inter in intermediates_LSU:
        LSU_up_dict[inter] = []; LSU_down_dict[inter] = []

        for index, rxn in assemblyRxns.iterrows():
            if rxn['product'] == inter:
                LSU_up_dict[inter].append(rxn['substrate'])

            if rxn['intermediate'] == inter:
                LSU_down_dict[inter].append(rxn['substrate'])
                 
                  

    return LSU_up_dict, LSU_down_dict

def getSSUSpecies(SSUfile_dir, N_inters):

    """
    Only return the set of species in SSU assembly pathways from 16S to SSU
    """
    def count_s(string):
        return string.count('s')
    
    smallAssemblyFile = SSUfile_dir
    assemblyData = json.load(open(smallAssemblyFile))
    mses = assemblyData['netmin_rmse']
    spsRemoved = assemblyData['netmin_species']

    spNames = set([ sp['name'] for sp in assemblyData['species'] if sp['name'][0] == 'R' ])

    maxImts = N_inters

    for err,sps in zip(mses,spsRemoved):
        if len(spNames) <= maxImts:
            break
        spNames.difference_update(set(sps))

    species_SSU = sorted(spNames, key=count_s)
    species_SSU[0] = '16S'

    
    return species_SSU




def getSSUInters(SSUfile_dir, N_inters):
    """
    Return: 
    substrates: Ordered List of substrates in assembly reactions containing N_inters intermediates
    intermediates_SSU: Ordered List of intermediates in assembly reactions containing N_inters intermediates
    products_SSU: Ordered List of products in assembly reactions containing N_inters intermediates
    """

    smallAssemblyFile = SSUfile_dir
    assemblyData = json.load(open(smallAssemblyFile))
    mses = assemblyData['netmin_rmse']
    spsRemoved = assemblyData['netmin_species']

    spNames = set([ sp['name'] for sp in assemblyData['species'] if sp['name'][0] == 'R' ])

    maxImts = N_inters

    for err,sps in zip(mses,spsRemoved):
        if len(spNames) <= maxImts:
            break
        spNames.difference_update(set(sps))

    assemblyRxns = [r for r in assemblyData['reactions'] if r['intermediate'] in spNames and r['product'] in spNames ]

    print(f"{len(assemblyRxns)} reactions connect {N_inters} intermediates in SSU")

    # print(assemblyRxns)

    # Sort substrates and intermediates consistently

    substrates = ['S' + rxn['protein'][1:] for rxn in assemblyRxns]

    intermediates_SSU = [rxn['intermediate'] for rxn in assemblyRxns]

    products_SSU = [rxn['product'] for rxn in assemblyRxns]
    
    metric = [_.count('s') for _ in products_SSU]

    sorted_indices = sorted(range(len(intermediates_SSU)), key=lambda i: metric[i])

    substrates = [substrates[i] for i in sorted_indices]

    intermediates_SSU = [intermediates_SSU[i] for i in sorted_indices]

    products_SSU = [products_SSU[i] for i in sorted_indices]

    # print(substrates, intermediates_SSU)

    # intermediates_SSU = sorted(spNames, key=count_s)
    # for a, b in zip(substrates, intermediates_SSU):
    #     print(a,b)
    
    intermediates_SSU[0] = '16S'

    return substrates, intermediates_SSU, products_SSU


def getLSUInters(LSUfile_dir):

    assemblyRxns = pd.read_excel(LSUfile_dir, sheet_name='LSU reactions', index_col=None)

    # for index, rxn in assemblyRxns.iterrows():
    #         intermediates_LSU.extend([rxn['intermediate'], rxn['product']])

    # intermediates_LSU = list(set(intermediates_LSU))

    # def count_L(string):
    #     return string.count('L') + string.count('S')

    # intermediates_LSU = sorted(intermediates_LSU, key = count_L)

    substartes_LSU = list(assemblyRxns['substrate'].values)
    intermediates_LSU = list(assemblyRxns['intermediate'].values)
    products_LSU = list(assemblyRxns['product'].values)

    intermediates_LSU[0] = '23S'
    
    return substartes_LSU, intermediates_LSU, products_LSU


def getCountInters(w, intermediates_list, SSULSUflag):
    """
    Return the count traces of intermediates of SSU or LSU
    """
    count_inters = np.zeros((len(intermediates_list), w.Nt, w.N_reps), dtype = np.float32)

    if SSULSUflag == 'SSU':
        for i_inter, inter in enumerate(intermediates_list):
                if inter == '16S':
                    count_R0069 = w.get_specie_trace('R_0069'); count_R0534 = w.get_specie_trace('R_0534')
                    count_16S = count_R0069 + count_R0534
                    count_inters[i_inter,:,:] = count_16S
                else:
                    count_inter = w.get_specie_trace(inter)
                    count_inters[i_inter,:,:] = count_inter

    elif SSULSUflag == 'LSU':
        for i_inter, inter in enumerate(intermediates_list):
                if inter == '23S':
                    count_R0068 = w.get_specie_trace('R_0068'); count_R0533 = w.get_specie_trace('R_0533')
                    count_23S = count_R0068 + count_R0533
                    count_inters[i_inter,:,:] = count_23S
                else:
                    count_inter = w.get_specie_trace(inter)
                    count_inters[i_inter,:,:] = count_inter
    else:
         print('Please input correct flag SSU or LSU')

    return count_inters

def getIDsToLocusNums(ribotoLocusTags):
    """
    
    Return: two dictionaries mapping IDs of precursors to their locusNums, e.g. 16S to 0069 and 0534

    """
    SSU_IDsToLocusNums = {'16S':['0069', '0534']}
    LSU_IDsToLocusNums = {'5S':['0067', '0532'], '23S':['0068', '0533']}

    for rptn_ID, LocusTags in ribotoLocusTags.items():
            if 'S' in rptn_ID and rptn_ID != 'S2' and rptn_ID != 'S21' and rptn_ID != 'S18': 
                    # S2 and S21 not included in the SSU assembly; S18 form dimer with S6
                    SSU_IDsToLocusNums[rptn_ID] = []
                    for locusTag in LocusTags:
                            locusNum = locusTag.split('_')[1]
                            SSU_IDsToLocusNums[rptn_ID].append(locusNum)
            
            elif 'L' in rptn_ID: # All 32 r-proteins are in assembly pathway
                    LSU_IDsToLocusNums[rptn_ID] = []
                    for locusTag in LocusTags:
                            locusNum = locusTag.split('_')[1]
                            LSU_IDsToLocusNums[rptn_ID].append(locusNum)

    return SSU_IDsToLocusNums, LSU_IDsToLocusNums


def getCountAvailable(w, IDsToLocusNums, SSULSUflag):

    count_availables = np.zeros((len(IDsToLocusNums), w.Nt, w.N_reps), dtype = np.float32)
    
    i_specie = 0
    
    if SSULSUflag == 'SSU':
            initialLetter = 'S'
    elif SSULSUflag == 'LSU':
            initialLetter = 'L'
    else:
            print('Please input correct flag SSU or LSU')
          
    for ID, LocusNums in IDsToLocusNums.items():
            temp = np.zeros((w.Nt, w.N_reps), dtype = np.float32)
            if ID.startswith(initialLetter): #r-proteins
                    for locusNum in LocusNums:
                            rptn  = 'P_' + locusNum
                            rptn_trace = w.get_specie_trace(rptn)
                            initial = rptn_trace[0,:] # 1D
                            Produced_rptn = 'Produced_P_' + locusNum
                            Produced_rptn_trace = w.get_specie_trace(Produced_rptn) #2D

                            temp = temp + initial + Produced_rptn_trace
            else: #rRNAs
                    for locusNum in LocusNums: 
                            rRNA  = 'R_' + locusNum
                            rRNA_trace = w.get_specie_trace(rRNA)
                            initial = rRNA_trace[0,:] # 1D
                            Produced_rRNA = 'Produced_R_' + locusNum
                            Produced_rRNA_trace = w.get_specie_trace(Produced_rRNA) #2D

                            temp = temp + initial + Produced_rRNA_trace


            count_availables[i_specie, :, :] = temp
            i_specie += 1

    return count_availables


def getBlockingIntervals(count_inters, inters_list, threshold):
    """
    
    Return:
    blockdict: intermediate name as key and value is the list of sublist. The sublist consists of replicate number and time intervals when above the threshold
    """

    block_dict = {}

    for i_inter, inter in enumerate(inters_list):
        count_inter = count_inters[i_inter, :, :]

        block_dict[inter] = []
        
        intervals = []
        
        for i_rep in range(count_inter.shape[1]):

            count_inter_rep = count_inter[:,i_rep]
            inter_above = count_inter_rep >= threshold

            if True in inter_above:
                intervals_rep = []                
                intervals_rep.append(i_rep)

                start = None
                for time_index, value_exceeds in enumerate(inter_above):
                    
                    if value_exceeds and start is None:
                        start = time_index
                    elif not value_exceeds and start is not None:
                        intervals_rep.append((start, time_index -1))
                        start = None
                if start is not None:
                    intervals_rep.append((start, len(inter_above)-1))
            
                intervals.append(intervals_rep)

            else:
                continue

        block_dict[inter] = intervals

    
    return block_dict


def getBlockInters(block_dict):
    
    blockInters = []

    for inter, intervals in block_dict.items():
        if intervals != []:
             blockInters.append(inter)

    return blockInters


def getBlockingInters_perReplicate(count_inters, inters_list, threshold):
    
    block_dict = {}

    for i_rep in range(count_inters.shape[2]):
        
        block_inters = []

        count_inters_rep = count_inters[:,:,i_rep]

        for i_inter, inter in enumerate(inters_list):
            
            count_inter_rep = count_inters_rep[i_inter, :]

            if np.max(count_inter_rep) > threshold:
                 
                block_inters.append(inter)
        
        block_dict[i_rep] = block_inters

    return block_dict

