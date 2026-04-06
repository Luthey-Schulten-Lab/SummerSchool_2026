"""

Description: functions to analyze ptns for CMEODE WCM
"""
import pandas as pd
import numpy as np
import WCM_gene as gene
import WCM_math as math

def gettRNAmap(genome):
    """
    Input: genome

    Description: 
    Set up the sub dictionary with amino acids and correspoding tRNAs {'LEU': ['R_0070', 'R_0423', 'R_0506'],...}
    """

    tRNA_map = {}

    for locusTag, locusDict in genome.items():
        
        if locusDict['Type'] == "tRNA":
            
            locusNum = locusTag.split('_')[1]
            
            rnaID = 'R_'+locusNum

            rnaName = locusDict['GeneName'].split('-')[1].upper()
            
            if rnaName not in tRNA_map:
                
                tRNA_map[rnaName] = [rnaID]
                
            else:
                
                tRNA_map[rnaName].append(rnaID)
                

    return tRNA_map


def getSynthetaseToAA(kineticfile_dir):
    synthetaseToAA = {}

    tRNA_params = pd.read_excel(kineticfile_dir, sheet_name='tRNA Charging')

    for index, row in tRNA_params.iterrows():
        if row['Parameter Type'] == 'synthetase':
            synthetase = row['Value']
            AA = row['Reaction Name'][0:3]
            synthetaseToAA[synthetase] = AA

    return synthetaseToAA


def gettRNALocusNums(tRNA_map):
    """
    
    Description: return a list of tRNA locusNums
    """
    tRNALocusNums = []


    return tRNALocusNums



def get_categories_ptn(init_conc_path):
    """
    Based on the proteomics in Excel
    Return: the list of locusNums that with experimental protemomic counts larger than 10 (non-ribosomal ptns)
        ptn_list contains 338 ptns, ribosomal_ptn_list contains 52, ten_list 65
    """


    ptn_list = []
    
    ribosomal_ptn_list = []
    
    ten_ptn_list = []

    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics')
    
    for row, protein in proteomics.iterrows():
        if row != 0:
            locusTag = protein['Locus Tag']
            locusNum = locusTag.split('_')[1]
            PtnName = protein['Gene Product']
    #         print(PtnName)
            if 'S ribosomal' not in PtnName: #exclude the ribosomal proteins
                initial_count = protein['Sim. Initial Ptn Cnt']
                if initial_count > 10:
                    ptn_list.append(locusNum)
                else:
                    ten_ptn_list.append(locusNum)
            else:
                ribosomal_ptn_list.append(locusNum)

    # print(f"ribosome Ptn: {ribosomal_ptn_list}")

    return ptn_list, ribosomal_ptn_list, ten_ptn_list

def getPtnLocations(init_conc_path):
    """
    
    Description: get the proteins at different locations
    """
    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics', skiprows=[1])

    locations = list(set(list(proteomics['Localization']))) 

    locations_ptns = {}

    for location in locations:
        locations_ptns[location] = []
        for index, row in proteomics.iterrows():
            if row['Localization'] == location:
                locusTag = row['Locus Tag']
                locusNum = locusTag.split('_')[1]
                locations_ptns[location].append(locusNum)

    for location in locations:
        print(f"{len(locations_ptns[location])} {location} proteins in Syn3A")

    sorted_dict = dict(sorted(locations_ptns.items(), key=lambda item: len(item[1]), reverse=True))

    return sorted_dict

def get_cplx_dict(cplx_path, locations_ptns):
    """
    Return a dictionary of complexes involved in the simulation

    Return: cplx_dict, two layer dictionary: {complex name, subdict}, {'init_count':0, 'Stoi': {'JCVISYN3A_0001':1, ...}, 'transMP_count': 2}
    """

    transmemPtnList = locations_ptns['trans-membrane']

    cplx_df = pd.read_excel(cplx_path, sheet_name='Complexes', dtype=str)
    # abc_df = pd.read_excel(cplx_path, sheet_name='ABC Transporter')
    # ATP_df = pd.read_excel(cplx_path, sheet_name='ATPSynthase', dtype=str)

    # cplx_dict = {}

    # def getABCTransporter(abc_complex, abc_df, subunit):

    #     locusNums = []
        
    #     for index, row in abc_df.iterrows():
    #         if row['Name'] == abc_complex:
    #             if row['Subunit'].startswith(subunit):
    #                 locusNums.append(row['Gene'])
        
    #     if locusNums != ['None']:
    #         locusTags = ['JCVISYN3A_' + locusNum for locusNum in locusNums]
    #     else:
    #         locusTags = []

    #     return locusTags



    # for index, row in cplx_df.iterrows():
    #     name = row['Name']
    #     cplx_dict[name] = {}
    #     # init_count
    #     cplx_dict[name]['init_count'] = row['Count']
    #     # Stoi
    #     type = row['Type']; cplx_dict[name]['Stoi'] = {}
    #     locusNums = row['Genes'].split(';'); locusTags = ['JCVISYN3A_' +locusNum for locusNum in locusNums]

    #     if type == 'dimer':
    #         for locusTag in locusTags:
    #             cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag)

    #     elif type == 'dimer_activate':
    #         for locusTag in locusTags:
    #             cplx_dict[name]['Stoi'][locusTag] = 1

    #     elif type == 'ECFModule':
    #         for locusTag in locusTags:
    #             cplx_dict[name]['Stoi'][locusTag] = 1

    #     elif type == 'ECF':
    #         cplx_dict[name]['Stoi'] = dict(cplx_dict['ECF']['Stoi'])
    #         cplx_dict[name]['Stoi'][locusTags[0]] = 1

    #     elif type == 'ABCTransporter':
    #         subunits = ['TMD', 'NBD', 'SBP']
    #         for subunit in subunits:
    #             locusTags = getABCTransporter(name, abc_df, subunit)
    #             for locusTag in set(locusTags):
    #                 cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag)

    #     elif type == 'ATPSynthase':
    #         for locusTag in locusTags:
    #             stoi = ATP_df[ATP_df['Gene'] == locusTag.split('_')[1]]['Stoi']
    #             cplx_dict[name]['Stoi'][locusTag] = int(stoi)
    #     elif type == 'RNAP':
    #         for locusTag in locusTags:
    #             cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag) 
    #         cplx_dict[name]['Stoi']['JCVISYN3A_0645'] = 2
    #     elif type == 'SecYEG':
    #         for locusTag in locusTags:
    #             cplx_dict[name]['Stoi'][locusTag] = locusTags.count(locusTag)
    #     else:
    #         print(f'WARNING: Unknown Assembly Pathway {type}')

    cplx_dict = {}


    for index, row in cplx_df.iterrows():
        name = row['Name']
        cplx_dict[name] = {}
        # init_count
        cplx_dict[name]['init_count'] = int(row['Init. Count'])
        # Stoi
        cplx_dict[name]['Stoi'] = {}

        locusNums = row['Genes Products'].split(';')
        Stois = [int(_) for _ in row['Stoichiometries'].split(';')]

        cplx_composition = []

        for locusNum, stoi in zip(locusNums, Stois): # get the composition of the cplx
            cplx_composition.extend([locusNum]*stoi)
        
        for locusNum in locusNums:
            cplx_dict[name]['Stoi'][locusNum] = cplx_composition.count(locusNum)

    for name, subdict in cplx_dict.items():
        Stoi = subdict['Stoi']
        transMP_count = 0
        for locusNum, stoi in Stoi.items():
            if locusNum in transmemPtnList:
                transMP_count += stoi
        subdict['transMP_count'] = transMP_count


    return cplx_dict



def getPtnInitCount(init_conc_path):
    """
    Get the initial count of ptn from Excel sheet
    """

    locusNumtoPtnInitCount = {}

    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics', skiprows=[1])

    for index, row in proteomics.iterrows():
        locusTag = row['Locus Tag']
        locusNum = locusTag.split('_')[1]
        init_count = row['Sim. Initial Ptn Cnt']
        locusNumtoPtnInitCount[locusNum] = init_count

    # print(f"{len(locusNumtoPtnInitCount)} Proteins in proteomics Excel Sheet")

    return locusNumtoPtnInitCount



def correctInitPtnCount(PtnIniCount, locusNum, cplx_dict, ribosomal_ptn_list, print_flag=False):
    """
    Input: 
        PtnIniCount: Count in proteomics Excel sheet
        locusTag
        cplx_dict: Python dictionary of complexes
        rPtn_locusTags

    Correct the initial protein counts for protein in complexes

    Return: 
        corrected_PtnIniCount: count of free protein after correction
        total_PtnIniCount: count of total ptn in the initial state after correction
        complex: complex that contains this protein 
    """


    complex = []

    total_PtnIniCount = 0
    
    old_PtnIniCount = PtnIniCount

    corrected_PtnIniCount = PtnIniCount


    # Assumption about the counts of free ribosomal protein pool (both SSU and LSU) is 5% of the complete assembled ribosome
    if locusNum in ribosomal_ptn_list:
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
        
        print(f"P_{locusNum}  Total Initial Count is {total_PtnIniCount}")

    return corrected_PtnIniCount, total_PtnIniCount, complex

def getPtnlist(init_conc_path, category):

    normal_ptn_list, ribosomal_ptn_list, ten_ptn_list = get_categories_ptn(init_conc_path)

    # rPtn_locusTags = ['JCVISYN3A_' + locusNum for locusNum in ribosomal_ptn_list]

    locusNumtoPtnInitCount = getPtnInitCount(init_conc_path)
    
    # ptns_lists = [normal_ptn_list, ribosomal_ptn_list, ten_ptn_list]

    if category == 'normal':
        ptns_list = normal_ptn_list
    elif category == 'rPtn':
        ptns_list = ribosomal_ptn_list
    elif category == 'ten':
        ptns_list = ten_ptn_list
    elif category == 'trans-membrane':
        TM_locusTags = getPtnLocations(init_conc_path)['trans-membrane']
        ptns_list = TM_locusTags
    elif category == 'entire':
        ptns_list = locusNumtoPtnInitCount.keys()

    return ptns_list
    

def getScaledPtnRatio(w, init_conc_path, locusNums, cplx_dict):
    """
    Input:
        init_conc_path: path to the proteomics Excel sheet
        locusNums: list of locusNums of proteins that fall into category
        ptns_ratio: 3D numpy array and the values are the ratio of total protein counts over each initial total protein count     
        surface_doubling_times: list of SA doubling time

    Output:

    Description:
        Check the scaled protein abundance at the SA doubling times
    """


    locusNumtoPtnInitCount = getPtnInitCount(init_conc_path)
    
    ribosomal_ptn_list = getPtnlist(init_conc_path, 'rPtn')

    Produced_ptns = ['Produced_P_{0}'.format(locusNum) for locusNum in locusNums ]

    Produced_ptns_counts = w.get_species_traces(Produced_ptns)

    ptns_ratio = np.zeros_like(Produced_ptns_counts)

    for i_ptn, locusNum in enumerate(locusNums):

        Produced_ptn_count = Produced_ptns_counts[i_ptn,:,:]

        PtnInitCount = locusNumtoPtnInitCount[locusNum]

        corrected_PtnIniCount, total_PtnIniCount, complex = correctInitPtnCount(PtnInitCount, locusNum, cplx_dict, ribosomal_ptn_list)

        scaled_produced_ptn_count = Produced_ptn_count / total_PtnIniCount

        scaled_total_ptn_count = scaled_produced_ptn_count + 1

        ptns_ratio[i_ptn,:,:] = scaled_total_ptn_count

    return ptns_ratio


def getAbnormalPtn(init_conc_path, locusNums, ptns_ratio, surface_doubling_times, cplx_dict, genomeDict, threshold=[1.5,3]):

    def getCplxcounts(complex):
        cplx_count = []

        for cplx in complex:
            count = cplx_dict[cplx]['init_count']
            cplx_count.append(count)

        return cplx_count

    def getStois(complex, locusNum):
        Stois = []

        for cplx in complex:
            Stoi = cplx_dict[cplx]['Stoi'][locusNum]
            Stois.append(Stoi)

        return Stois
        
    abnormal_ptns = []; abnormal_ratios = []
    
    abnormal_ptns_dicts = []

    # get rPtn_locusTags
    ribosomal_ptn_list = getPtnlist(init_conc_path, 'rPtn')

    proteomics = pd.read_excel(init_conc_path, sheet_name='Comparative Proteomics', skiprows=[1])

    for i_ptn, locusNum in enumerate(locusNums):
        ptn_ratio = ptns_ratio[i_ptn, :, :] # 2D array

        min, median, mean, max = math.get_min_median_mean_max(math.get_doubling_moments_value(surface_doubling_times,ptn_ratio)) # min, median, mean, max of scaled protein counts

        if mean < threshold[0] or mean > threshold[1]:
            abnormal_ptns.append(locusNum); abnormal_ratios.append(mean)

            abnormal_ptn_dict = {}
            locusTag = 'JCVISYN3A_' + locusNum
            row = proteomics[proteomics['Locus Tag'] == locusTag]
            PtnIniCount = row['Sim. Initial Ptn Cnt'].values[0]

            corrected_PtnIniCount, total_PtnIniCount, complex = correctInitPtnCount(PtnIniCount, locusNum, cplx_dict, ribosomal_ptn_list, print_flag=False)
            
            abnormal_ptn_dict['Locus Tag'] = locusTag
            abnormal_ptn_dict['Gene Name'] = row['Gene Name'].values[0]
            abnormal_ptn_dict['Gene Product'] = row['Gene Product'].values[0]
            abnormal_ptn_dict['Protein Length'] = len(genomeDict[locusTag]['AAsequence'])
            abnormal_ptn_dict['Localization'] = row['Localization'].values[0]

            abnormal_ptn_dict['Complex'] = ', '.join(a for a in complex)
            cplx_count = getCplxcounts(complex)
            abnormal_ptn_dict['Sim. Initial Cplx Cnt'] = ', '.join(str(a) for a in cplx_count)
            Stois = getStois(complex, locusNum)
            abnormal_ptn_dict['Stois in Cplx'] = ', '.join(str(a) for a in Stois)

            abnormal_ptn_dict['Exp. Ptn Cnt'] = row['Exp. Ptn Cnt'].values[0]

            # ptn['Sim. Initial Ptn Cnt'] = PtnIniCount
            # abnormal_ptn_dict['Sim. Free Initial Ptn Cnt'] = corrected_PtnIniCount
            # abnormal_ptn_dict['Sim. Total Initial Ptn Cnt'] = total_PtnIniCount
            
            abnormal_ptn_dict['w/o. Cplx Ptn Cnt'] = row['Sim. Initial Ptn Cnt'].values[0]
            abnormal_ptn_dict['w. Cplx Free Ptn Cnt'] = corrected_PtnIniCount
            abnormal_ptn_dict['w. Cplx Total Ptn Cnt'] = total_PtnIniCount

            abnormal_ptn_dict['Min Scaled Abund'] = f"{min:.2f}"
            abnormal_ptn_dict['Median Scaled Abund'] = f"{median:.2f}"
            abnormal_ptn_dict['Max Scaled Abund'] = f"{max:.2f}"
            abnormal_ptn_dict['Mean Scaled Abund'] = f"{mean:.2f}"

            abnormal_ptns_dicts.append(abnormal_ptn_dict)
    
    abnornal_ptn_df = pd.DataFrame(abnormal_ptns_dicts)

    # # print(abnornal_ptn_df)
    # if abnornal_ptn_df.empty:
    #     sorted_abnornal_ptn_df = pd.DataFrame()
    # else:
    #     sorted_abnornal_ptn_df = abnornal_ptn_df.sort_values(by='Mean Scaled Abund')

    return abnormal_ptns, abnormal_ratios, abnornal_ptn_df
    
def getRxnsDict():
    """
    Description: for the old kinetic_params Excel sheet
    {'BPNT': {'locusNums': ['0139']},
    'FAKr': {'locusNums': ['0420', '0616', '0617'], 'GPR': 'and'},...}

    """
    sheet_names = ['Central', 'Nucleotide', 'Lipid', 
                    'Cofactor', 'Transport']
    kinetic_params = '/data/enguang/CMEODE/input_data/kinetic_params.xlsx'
    rxns_dict = {}
    for sheet_name in sheet_names:
        params = pd.read_excel(kinetic_params, sheet_name=sheet_name)
        for index, row in params.iterrows():
            if row['Parameter Type'] == 'Eff Enzyme Count':
                rxn_id = row['Reaction Name']
                ptns = row['Value']
                if ptns != 'default':
                    ptns = ptns.split('-')
                    locusNums = [ptn.split('_')[1] for ptn in ptns]
                    rxns_dict[rxn_id] = {}
                    rxns_dict[rxn_id]['GPR_locusNums'] = locusNums
                    
                    if len(locusNums) > 1:
                        rxn_params = params[params['Reaction Name'] == rxn_id]
                        GPRrule = rxn_params[rxn_params['Parameter Type'] == 'GPR rule']['Value'].values[0]
                        rxns_dict[rxn_id]['GPR'] = GPRrule

                        
                else:
                    print(f"Reaction {rxn_id} has no assigned proteins")

    return rxns_dict

def getChangedRxnsDict(rxns_dict, cplx_dict):
    """
    
    Get the information of metabolic reactions that changed from GPR AND/OR rule to complex
    """
    changed_rxns_dict = {}

    print("********* Checking Changed Reactions ***********")

    for rxn_id, rxn_subdict in rxns_dict.items():
        GPR_locusNums = rxn_subdict['GPR_locusNums']

        for cplx_id, cplx_subdict in cplx_dict.items():
            cplx_stoi = cplx_subdict['Stoi']
            cplx_locusTags = cplx_stoi.keys()
            cplx_locusNums = [locusTag.split('_')[1] for locusTag in cplx_locusTags]
            if cplx_id != 'ECF':
                if set(GPR_locusNums).issubset(set(cplx_locusNums)):
                    print(f"{rxn_id} Old reaction locusNums {GPR_locusNums} cplx {cplx_id} {cplx_locusNums}")
                    changed_rxns_dict[rxn_id] = rxn_subdict
                    changed_rxns_dict[rxn_id]['complex'] = cplx_id
                    
                elif set(cplx_locusNums).issubset(set(GPR_locusNums)):
                    print(f"{rxn_id} Old reaction locusNums {GPR_locusNums} cplx {cplx_id} {cplx_locusNums}")
                    changed_rxns_dict[rxn_id] = rxn_subdict
                    changed_rxns_dict[rxn_id]['complex'] = cplx_id

    print(changed_rxns_dict)

    return changed_rxns_dict


def checkComplexStoi(w, surface_doubling_times, cplx_dict, cplx, produced_prefix = 'Produced_P_'):
    """
    
    Return:
        ratio_dict: {'0009': 1, '0010': 2}

    Description: 
        Check the Stoichiometric balance among generated subunits up to SA doubling time; 
        Divide the generated subunits counts over the initial complex count
    """

    ratio_dict = {}

    Stoi = cplx_dict[cplx]['Stoi']

    init_cplx_count = cplx_dict[cplx]['init_count'] # initial count of assembled complex

    print(f"Stoi of complex {cplx}: {Stoi}")

    locusNums = list(Stoi.keys())

    produced_ptns_ids = [produced_prefix + locusNum for locusNum in locusNums ]

    produced_subunits_counts = w.get_species_traces(produced_ptns_ids)

    counts_SA = math.get_doubling_moments_value(surface_doubling_times, produced_subunits_counts)

    ratio_SA = counts_SA/init_cplx_count

    avg_ratio_SA = np.mean(ratio_SA, axis=1)

    for locusNum, ratio in zip(locusNums, avg_ratio_SA):
        stoi = Stoi[locusNum]
        ratio_dict[locusNum] = f"{ratio/stoi:.2f}"

    ratio_dict['init_count'] = init_cplx_count

    return ratio_dict
 

def plotProducedPtns(w, surface_doubling_times, cplx_dict, cplx, fig_dir, fig_label, reps, produced_prefix = 'Produced_P_'):
    """

    Description: plot the traces and histogram of produced protein subunits of a complex 
    """

    Stoi = cplx_dict[cplx]['Stoi']

    init_cplx_count = cplx_dict[cplx]['init_count'] # initial count of assembled complex

    print(f"Stoi of complex {cplx}: {Stoi}")

    locusTags = list(Stoi.keys())

    produced_ptns_ids = [produced_prefix + locusTag.split('_')[1] for locusTag in locusTags ]

    produced_subunits_counts = w.get_species_traces(produced_ptns_ids)

    counts_SA = math.get_doubling_moments_value(surface_doubling_times, produced_subunits_counts)

    for i_subunit, subunit in enumerate(locusTags):

        ylabel = 'Count [\#]'
        title = produced_ptns_ids[i_subunit] + f' in {cplx} with init count {init_cplx_count} and Stoi {Stoi[subunit]}'
        w.plot_in_replicates_single(fig_dir,fig_label,
                               '.png',
                              produced_subunits_counts[i_subunit,:,:], reps, ylabel, title,
                               True, True)
        
        xlabel = produced_ptns_ids[i_subunit] + ' at SA doubling time'
        produced_subunits_SA = counts_SA[i_subunit,:]
        w.plot_hist(fig_dir, fig_label, '.png', produced_subunits_SA, xlabel, ylabel, title, bins=50)

    return None

def getSubcomplexes(w, cplx_dict):
    """
    Return:
        sub_cplx_lists
    Description: 
        get the lists of subcomplexes under one holo-complex from the names of species
    """

    species = [_.replace('Produced_','') for _ in w.species]
    species = list(set(species))

    sub_cplx_lists = {}
    for holo_cmplx in cplx_dict.keys():
        sub_cplx_lists[holo_cmplx] = []

    for specie in species:
        for holo_cmplx in cplx_dict.keys():
            if specie.startswith(holo_cmplx+'_'):
                if holo_cmplx == 'Degradosome':
                    if specie.find('mRNA_') == -1:
                        sub_cplx_lists[holo_cmplx].append(specie)
                elif holo_cmplx == 'RNAP':
                    if specie.find('G_') == -1:
                        sub_cplx_lists[holo_cmplx].append(specie)
                elif holo_cmplx == 'SecYEGDF':
                    if specie.find('RNC_') == -1:
                        sub_cplx_lists[holo_cmplx].append(specie)
                else:
                    sub_cplx_lists[holo_cmplx].append(specie)

    sub_cplx_lists['SecYEGDF'] = ['P_0652', 'SecYE', 'SecYEG']

    # Sort the strings by length
    for complex, intermediates in sub_cplx_lists.items():
        sorted_intermediates = sorted(intermediates, key=len)
        sub_cplx_lists[complex] = sorted_intermediates

    return sub_cplx_lists

def getCplxIntermediatesMap(cplx_dict, cplx_path, sub_cplx_lists, transmemPtnList):
    """
    Description: 
        Get the protein complex intermediates/subcomplex's dictionary

    Assume initial counts of complex intermediates are all 0
    """

    holo_cmplx = list(cplx_dict.keys())
    produced_holo_cmplx = ['Produced_'+_ for _ in holo_cmplx]

    letter_locusNum = getLettertolocusNum()

    sub_cplx_dict = {}
    
    for cplx in sum(sub_cplx_lists.values(), []):
        if cplx not in holo_cmplx+produced_holo_cmplx: # Exclude the completely assembled complex
            sub_cplx_dict[cplx] = {}

            sub_cplx_dict[cplx]['init_count'] = 0 # Always 0

            Stoi = mapStringtoStoi(cplx, cplx_path, letter_locusNum) # 

            sub_cplx_dict[cplx]['Stoi'] = Stoi
    

    for name, subdict in sub_cplx_dict.items():
        Stoi = subdict['Stoi']
        transMP_count = 0
        for locusNum, stoi in Stoi.items():
            if locusNum in transmemPtnList:
                transMP_count += stoi
        subdict['transMP_count'] = transMP_count
        
    # print(f"sub_cplx_dict {len(sub_cplx_dict)}")

    # for cplx, subdict in sub_cplx_dict.items():
    #     print(cplx, subdict)

    return sub_cplx_dict

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


def getLocusNumRNCs(locations_ptns, YidC_ptns):
    """
    Description:
        Get the names of species of RNCs of a certain protein locusNum
    """
    import itertools

    LocusNum2RNCs = {}
    
    locusNums = list(itertools.chain(*locations_ptns.values()))

    for locusNum in locusNums:
        if locusNum in locations_ptns['peripheral membrane']+locations_ptns['cytoplasm']+locations_ptns['unidentified']:
            RNC_prefixes = ['RNC_']
        elif locusNum in locations_ptns['extracellular']:
            RNC_prefixes = ['RNC_', 'SecA_RNC_', 'SecYEGDF_SecA_RNC_']
        elif locusNum in locations_ptns['lipoprotein']:
            RNC_prefixes = ['RNC_', 'SecA_RNC_', 'SecYEGDF_SecA_RNC_', 'RNC_Long_']
        elif locusNum in locations_ptns['trans-membrane']:
            if locusNum in YidC_ptns:
                RNC_prefixes = ['RNC_']
            else:
                RNC_prefixes = ['RNC_', 'SRP_RNC_', 'SR_SRP_RNC_', 'SecYEGDF_RNC_', 'RNC_Long_']

        RNCs = [_ + locusNum for _ in RNC_prefixes]
        LocusNum2RNCs[locusNum] = RNCs

    return LocusNum2RNCs


def getLocation2RNCs(locations_ptns, LocusNums2RNCs, YidC_ptns):
    """
    Description:
        Get the names of RNCs in certain locations
    """
    
    Location2RNCs = {}
    Location2RNCs['Cytoplasmic']  = []
    Location2RNCs['SRP']  = []
    Location2RNCs['SecA']  = []

    cyto_LocusNums = locations_ptns['peripheral membrane'] + locations_ptns['cytoplasm']+locations_ptns['unidentified']+YidC_ptns
    
    SecA_LocusNums = locations_ptns['secreted'] + locations_ptns['lipoprotein']

    SRP_LocusNums = list(set(locations_ptns['trans-membrane']) - set(YidC_ptns))

    for locusNum in cyto_LocusNums:

        Location2RNCs['Cytoplasmic'].extend(LocusNums2RNCs[locusNum])

    for locusNum in SRP_LocusNums:
        Location2RNCs['SRP'].extend(LocusNums2RNCs[locusNum])

    for locusNum in SecA_LocusNums:
        Location2RNCs['SecA'].extend(LocusNums2RNCs[locusNum])

    return Location2RNCs

def getPathway2Secs(w, SRP_LocusNums, SecA_LocusNums):
    """
    
    Description:
        Get the species in the SRP and SecA pathways
    """

    Pathway2Secs = {}
    Pathway2Secs['SRP'] = []; Pathway2Secs['SecA'] = []

    for specie in w.species:
        if specie.startswith('SecYEGDF_'):
            if specie[-4:] in SRP_LocusNums:
                Pathway2Secs['SRP'].append(specie)
            elif specie[-4:] in SecA_LocusNums:
                Pathway2Secs['SecA'].append(specie)
            else:
                print(specie)

    return Pathway2Secs

def cal_translation_per_mRNA(w, locusNums, reps, surface_doubling_times, locations_ptns, lipo=True):
    """
    Input:
        lipo:
            True: use Produced_P_
            False: use Produced_PreP
    Description: 
        Calculate the translation events per mRNA
    """
    
    trans_per_mRNA = []
    
    zero_mRNA_list = []
    # mRNA_prefix = 'R_'
    produced_mRNA = 'Produced_R_'
    produced_ptn = 'Produced_P_'
    produced_pre_ptn = 'Produced_PreP_'
    lipo_locusNums = locations_ptns['lipoprotein']

    for locusNum in locusNums:
        if locusNum in lipo_locusNums and not lipo: # use Produced_PreP_
            produced_prefix = produced_pre_ptn
        else:
            produced_prefix = produced_ptn

        trans_reps = []
        for rep, SA_time in zip(reps, surface_doubling_times):
            
            P_R = w.get_specie_trace(produced_mRNA+locusNum)[int(SA_time),rep-1]

            P_P = w.get_specie_trace(produced_prefix+locusNum)[int(SA_time),rep-1]
            
            if P_R != 0:
                trans_reps.append(P_P/P_R)
            else:
                zero_mRNA_list.append(locusNum)
        

        trans_per_mRNA.append(np.mean(trans_reps))

    return  np.array(trans_per_mRNA), list(set(zero_mRNA_list))

def get_produced_subunits_holo_cplx(w, TypetoLocusNums, stois, locusNums, cplx):
    """
    Input:
        TypetoLocusNums:
            Gene types to locusNums
    Description:
        Get the produced counts of protein subunits in cplx divied by their stoichiometry
    """
    
    ptn_locusNums = TypetoLocusNums['protein']

    ptn_prefix = 'Produced_P_'; rRNA_tRNA_prefix = 'Produced_R_'

    produced_species = [ptn_prefix + _ if _ in ptn_locusNums else rRNA_tRNA_prefix + _ for _ in locusNums]

    produced_species.append('Produced_'+cplx)

    produced_counts = math.get_doubling_moments_value(w.surface_doubling_times, w.get_species_traces(produced_species)) # 2D Specie by replicates

    for i, locusNum in enumerate(locusNums):
        produced_counts[i,:] = produced_counts[i,:]/stois[locusNum]

    return produced_counts

def get_ava_subunits(w, TypetoLocusNums, stois, locusNums, cplx):
    """
    
    Description:
        Return the available count of subunits without produced complex
    """

    ptn_locusNums = TypetoLocusNums['protein']

    ptn_prefix = 'Produced_P_'; rRNA_tRNA_prefix = 'Produced_R_'

    produced_species = [ptn_prefix + _ if _ in ptn_locusNums else rRNA_tRNA_prefix + _ for _ in locusNums]

    produced_counts = w.get_species_traces(produced_species) # 3D Specie by time by replicates (N,T,reps)
    init_subunits_count = w.get_species_traces(['P_'+_ for _ in locusNums])[:,0:1,0:1] # (N,1,1) array

    ava_counts = produced_counts + init_subunits_count

    for i, locusNum in enumerate(locusNums):
        ava_counts[i,:] = ava_counts[i,:]/stois[locusNum] # Over Stoi


    return ava_counts # 3D array

def get_bottleneck_species(w, TypetoLocusNums, stois, locusNums, cplx):
    """
    For non-ribosome complexes since 23S, 5S and 16S rRNAs are encoded by two genes
    """
    ava_counts = get_ava_subunits(w, TypetoLocusNums, stois, locusNums, cplx) # 3D array

    ava_counts_SA = math.get_doubling_moments_value(w.surface_doubling_times,
                                                    ava_counts)

    subuits_indices = np.argmin(ava_counts_SA, axis = 0) # length of # of replicates

    bottle_neck_locusNums = [locusNums[index] for index in subuits_indices]

    gene_names = gene.get_gene_names('/data/enguang/CMEODE/input_data/initial_concentrations_AfterDavid.xlsx', bottle_neck_locusNums, w.genomeDict)

    bottle_neck_locus_gene = [a+'/'+b for a,b in zip(gene_names, bottle_neck_locusNums)]    

    return bottle_neck_locusNums, bottle_neck_locus_gene

# def plotIntermediatesCplx(w, ):


#     return None


# def compareCplx(cplx_dict, rxns_dict):

#     compared_cplx_dict = {}

    
    
#     return compared_cplx_dict