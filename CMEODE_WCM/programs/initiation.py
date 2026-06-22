"""
Author: Enguang Fu

Date: March 2024

initialize simulation conditions, necessary constants 
"""


import pandas as pd
import numpy as np
from Bio.Seq import Seq
from itertools import chain


import rxns_CME

#########################################################################################
# Mapping gene bank file into Dict
def mapDNA(genome):
    DNAmap = {}
    chromosome_length = int((genome.features[0].location.end+1)/10)
    ori_ter_rotation_factor = int(chromosome_length/2)
    print(ori_ter_rotation_factor)   
    for feature in genome.features:
        strand = feature.strand     
        if strand == 1:    
            start = int(feature.location.start/10)        
            if start<=ori_ter_rotation_factor:
                start = start + ori_ter_rotation_factor           
            elif start>ori_ter_rotation_factor:               
                start = start - ori_ter_rotation_factor
            end = int(feature.location.end/10)
            if end<=ori_ter_rotation_factor:             
                end = end + ori_ter_rotation_factor
            elif end>ori_ter_rotation_factor:
                end = end - ori_ter_rotation_factor
        elif strand == -1:      
            start = int(feature.location.end/10)        
            if start<=ori_ter_rotation_factor:
                start = start + ori_ter_rotation_factor
            elif start>ori_ter_rotation_factor:
                start = start - ori_ter_rotation_factor
            end = int(feature.location.start/10)
            if end<=ori_ter_rotation_factor:            
                end = end + ori_ter_rotation_factor        
            elif end>ori_ter_rotation_factor:
                end = end - ori_ter_rotation_factor
        if feature.type == 'CDS':
            
            if('protein_id' in feature.qualifiers.keys()):
                # Excluding 3 pseudo genes
                # Totoal coding genes are 455
                locusTag = feature.qualifiers['locus_tag'][0]
                DNAmap[locusTag] = {}
                DNAmap[locusTag]['Type'] = 'protein'
                DNAmap[locusTag]['startIndex'] = [int(start)]
                DNAmap[locusTag]['originalStart'] = int(start)
                DNAmap[locusTag]['endIndex'] = [int(end)]           
                DNAmap[locusTag]['originalEnd'] = int(end)
                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
                DNAmap[locusTag]['AAsequence'] = str(feature.location.extract(genome.seq).transcribe().translate(table=4))            
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])     
        elif feature.type == 'tRNA':
            locusTag = feature.qualifiers['locus_tag'][0]       
            DNAmap[locusTag] = {} 
            DNAmap[locusTag]['Type'] = 'tRNA'
            DNAmap[locusTag]['startIndex'] = [int(start)]          
            DNAmap[locusTag]['originalStart'] = int(start)
            DNAmap[locusTag]['endIndex'] = [int(end)]
            DNAmap[locusTag]['originalEnd'] = int(end) 
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])      
        elif feature.type == 'rRNA':           
            locusTag = feature.qualifiers['locus_tag'][0]           
            DNAmap[locusTag] = {}          
            DNAmap[locusTag]['Type'] = 'rRNA'
            DNAmap[locusTag]['startIndex'] = [int(start)]
            DNAmap[locusTag]['originalStart'] = int(start)
            DNAmap[locusTag]['endIndex'] = [int(end)]
            DNAmap[locusTag]['originalEnd'] = int(end)          
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            
        else:
            # 2 misc binding features, 2 ncRNA genes, 1 tm RNA gene, 
            try:
                locusTag = feature.qualifiers['locus_tag'][0]           
                DNAmap[locusTag] = {}          
                DNAmap[locusTag]['Type'] = feature.type
                DNAmap[locusTag]['startIndex'] = [int(start)]
                DNAmap[locusTag]['originalStart'] = int(start)
                DNAmap[locusTag]['endIndex'] = [int(end)]
                DNAmap[locusTag]['originalEnd'] = int(end)          
                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            except:
                continue
            
    return DNAmap, chromosome_length
#########################################################################################

def categorizeGenes(sim_properties):

    
    genomeDict = sim_properties['genome']

    geneTypes = []
    geneNumbers = []
    TypetolocusNums = {}

    for locusTag, locusDict in genomeDict.items():
        locusNum = locusTag.split('_')[1]

        if locusDict['Type'] not in geneTypes:
            geneTypes.append(locusDict['Type'])
            TypetolocusNums[locusDict['Type']] = []

        TypetolocusNums[locusDict['Type']].append(locusNum)
    
    for type, locusNums in TypetolocusNums.items():
        geneNumbers.append(len(locusNums))

    # map(str, list) will convert the list's elements into strings.
    print('Six types of genes in Syn3A are {0} with respective numbers {1}.'.format(', '.join(geneTypes), ', '.join(map(str, geneNumbers))))
    # 'Type' = ['protein', 'ncRNA', 'gene', 'rRNA', 'tRNA', 'tmRNA']
    locusNumtoType = {}

    for type, locusNums in TypetolocusNums.items():
        for locusNum in locusNums:
            locusNumtoType[locusNum] = type
    
    sim_properties['TypetolocusNums'] = TypetolocusNums
    sim_properties['locusNumtoType'] = locusNumtoType

    return TypetolocusNums, locusNumtoType, geneTypes

#########################################################################################
def initializeCME(sim, restartNum, sim_properties):
    """
    
    Decription: Set up the simulation time and hook and write intervals for CME simulation
    """
    
    # Initialize every restart CME simulation object
    restartInterval = sim_properties['restartInterval']
    hookInterval = sim_properties['hookInterval']
   
    sim.setWriteInterval(hookInterval)
    sim.setHookInterval(hookInterval)
    sim.setSimulationTime(restartInterval)
    
    return sim
#########################################################################################



#########################################################################################
def initializeConstants(sim_properties):
    """
    
    Description: create constants including lists or sub dictionaries in sim_properties. 
        Sub dictionaries: tRNAMap, promoter strengths, and locusNumtoIndex
        Lists: membrane proteins list, pseodogenes list

    """

    kinetic_params_path =  '../input_data/kinetic_params.xlsx'; sim_properties['kinetic_params_path'] = kinetic_params_path
    init_conc_path = '../input_data/initial_concentrations.xlsx'; sim_properties['init_conc_path'] = init_conc_path
    cplx_path = '../input_data/complex_formation.xlsx'; sim_properties['cplx_path'] = cplx_path
    
    categorizeGenes(sim_properties)
    gettRNAMap(sim_properties)
    getPtnLocations(sim_properties); getPtnInitCount(sim_properties); getAbnormalTranslocationPtns(sim_properties)
    getComplexMap(sim_properties); getTransloconRNCComplex(sim_properties)
    getriboTolocusMap(sim_properties)
    initializePromoterStrengths(sim_properties)


    sim_properties['aa_list'] = ['M_ala__L_c','M_arg__L_c', 'M_asn__L_c', 'M_asp__L_c', 'M_cys__L_c',
        'M_glu__L_c', 'M_gly_c', 'M_his__L_c', 'M_ile__L_c', 'M_leu__L_c',
        'M_lys__L_c', 'M_met__L_c', 'M_phe__L_c', 'M_pro__L_c', 'M_ser__L_c',
        'M_thr__L_c', 'M_trp__L_c', 'M_tyr__L_c', 'M_val__L_c', 'M_gln__L_c']
    
    sim_properties['rnap_spacing'] = int(400)

    sim_properties['pseudoGenes'] = ['JCVISYN3A_0051', 'JCVISYN3A_0546','JCVISYN3A_0602']

    r_cell = 2e-7 # m
    sim_properties['r_cell'] = r_cell
    
    cellVolume_init = (4*np.pi/3)*r_cell**3*1000 # L
    sim_properties['volume_L'] = [cellVolume_init]

    print(f"The initial radius of syn3A is {r_cell*1e9:.2f} nm")

    # sim_properties['counts'] is a dictionary that record the number of all species per hookInterval
    sim_properties['counts'] = {}

    # sim_properties['conc'] is a dictionary that record the concentrations of metabolites after each ODE run
    # sim_properties['conc'] = {}

    # sim_properties['fluxes'] is a dictionary that record the fluxes of rxns after each ODE run
    sim_properties['fluxes'] = {}

    return None
#########################################################################################

#########################################################################################
def getPtnLocations(sim_properties):
    """
    
    Description: get the proteins at different locations
    """
    proteomics = pd.read_excel(sim_properties['init_conc_path'], sheet_name='Comparative Proteomics', skiprows=[1])

    locations = list(set(list(proteomics['Localization']))) 

    locations_ptns = {}

    for location in locations:
        locations_ptns[location] = []
        for index, row in proteomics.iterrows():
            if row['Localization'] == location:
                locusTag = row['Locus Tag']
                locusNum = locusTag.split('_')[1]
                locations_ptns[location].append(locusNum)
    
    sim_properties['locations_ptns'] = locations_ptns

    num_loca_ptns = [len(locusTags) for location, locusTags in locations_ptns.items()]

    # Sort the location string based on the number of proteins for better printing
    sorted_pairs = sorted(zip(locations, num_loca_ptns), key=lambda x: x[1], reverse=True)
    locations = [location for location, _ in sorted_pairs]
    num_loca_ptns = [num for _, num in sorted_pairs]

    print(f"{sum(num_loca_ptns)} proteins' localizatoin are categorized into {', '.join(locations)} with respective number {', '.join(map(str, num_loca_ptns))}")


    return None
#########################################################################################

#########################################################################################
def getPtnInitCount(sim_properties):
    """
    Get the initial count of ptn from Excel sheet
    """

    locusTagtoPtnInitCount = {}

    proteomics = pd.read_excel(sim_properties['init_conc_path'], sheet_name='Comparative Proteomics', skiprows=[1])

    for index, row in proteomics.iterrows():
        locusTag = row['Locus Tag']
        init_count = row['Sim. Initial Ptn Cnt']
        locusTagtoPtnInitCount[locusTag] = init_count
    
    sim_properties['proteomics_count'] = locusTagtoPtnInitCount
        
    return None
#########################################################################################

#########################################################################################
def getAbnormalTranslocationPtns(sim_properties):
    """
    Get the initial count of ptn from Excel sheet
    """

    df = pd.read_excel(sim_properties['cplx_path'], sheet_name='Sole YidC IMP', dtype=str)

    sim_properties['sole_YidC_ptns'] = list(df['locusNum'])

    df = pd.read_excel(sim_properties['cplx_path'], sheet_name='Peri MP in Cyto', dtype=str)

    sim_properties['peri_MP_in_cyto'] = list(df['locusNum'])
        
    return None
#########################################################################################


#########################################################################################
def initializePromoterStrengths(sim_properties):
    """
    Input: sim_properties

    Return: None

    Called at the beginning of the whole simulation once

    Description: 
    """

    # Promotoer Strength maintain in the whole cell cycle
    # Used to modify the transcription rates of mRNA, tRNA, and rRNA
    # Ribosomal proteins have higher transcription rate

    # kcat_mod: tRNA 85, rRNAs: 85
    PromoterStrengths = {}
    
    genome = sim_properties['genome']   

    proteomics = pd.read_excel(sim_properties['init_conc_path'], sheet_name='Comparative Proteomics', skiprows=[1])

    cplx_dict = sim_properties['cplx_dict']

    ATPSynthase_stoi = {'0789':1,'0790':3, '0791':1 , '0792':3, '0793':1, '0794':2, '0795':10, '0796':1}

    for locusTag, locusDict in genome.items():    
        locusNum = locusTag.split('_')[1]

        if locusDict["Type"] == 'protein':
            PtnInitCount = proteomics.loc[ proteomics['Locus Tag'] == locusTag ]['Sim. Initial Ptn Cnt'].values[0]

            corrected_PtnInitCount, total_PtnInitCount, complex = rxns_CME.correctInitPtnCount(PtnInitCount, locusNum, sim_properties)

            PromoterStrengths[locusNum] = (min(765, max(10, total_PtnInitCount)))/180
                         
        elif locusDict["Type"] == 'tRNA':
            # PromoterStrengths[locusNum] = 180
            # Different from cell paper
            PromoterStrengths[locusNum] = 765/180

        elif locusDict["Type"] == 'rRNA':
            # Different from cell paper
            # Since 23S is 2900 nt long so give a higher promoter strength

            if locusNum in ['0068', '0533']: # 23S
                PromoterStrengths[locusNum] = 765*9/180
            elif locusNum in ['0069', '0534']: # 16S
                PromoterStrengths[locusNum] = 765*6/180
            elif locusNum in ['0067', '0532']: # 5S
                PromoterStrengths[locusNum] = 765*2/180

        else:
            PromoterStrengths[locusNum] = 10/180
            
    sim_properties['promoters'] = PromoterStrengths

    print(f"Promoter Strength for transcription Initialized")

    return None
#########################################################################################

#########################################################################################
def initializeFluxes(sim_properties, odemodel, boolean_end):
    """
    Called after the first time odecell object built

    Description: create the flux dictionary to store the flux trajectories 
    """
    
    rxns = odemodel.getRxnList()

    sim_properties['fluxes'] = {}

    fluxesDict = sim_properties['fluxes']

    for rxn in rxns:
        rxn_id = rxn.getID()
        fluxesDict['F_' + rxn_id] = []
        if boolean_end:
            fluxesDict['F_' + rxn_id + '_end'] = []

    return None
#########################################################################################

#########################################################################################
# def initializeConcDictionary(sim_properties, odemodel):
#     """"
#     Called after the first time ODE run

#     Description: initialize the trajectories of metabolites in ODE 
#     """
#     concDict = sim_properties['conc']
#     metIDDict = odemodel.getMetDict()
#     # {'M_ACP_c': 0, 'M_ACP_R_c': 1, 'M_apoACP_c': 2, ...}

#     for species in metIDDict.keys():
#         concDict[species] = []

#     return None

#########################################################################################

#########################################################################################
def getReactionMap(sim_properties):
    """
    
    Description: initialize reaction map where the keys are the names and values are serial number, i.e. 1, 2, ...
    """
    
    sim_properties['rxns_map'] = {} # rxns_map {'init_1';1, 'init_2':2,...}

    rxns_prefix = {'initiation':'init', 'replication':'rep', 'trsc_binding':'trsc_binding', 'trsc_poly':'trsc_poly',
                   'tran_binding':'tran_binding', 'tran_poly':'tran_poly', 'translocation':'translocation','deg_binding':'deg_binding', 'deg_depoly': 'deg_depoly', 
                   'ribo_biogenesis':'ribo', 'tRNACharging':'tRNA', 'ptnComplex':'cplx'}  #Prefix for reaction in subsystems
    sim_properties['rxns_prefix'] = rxns_prefix

    rxns_numbers = {}
    for key in rxns_prefix.keys():
        rxns_numbers[key] = [int(0)]

    # rxns_numbers = {'initiation': [int(0)], 'replication':[int(0)], 'trsc_binding':[int(0)], 'trsc_ploy':[int(0)],
    #                'tran_binding':[int(0)], 'tran_poly':[int(0)], 'translocation':[int(0)], 'deg_binding':[int(0)], 'deg_depoly':[int(0)], 
    #                'ribo_biogenesis':[int(0)], 'tRNACharging':[int(0)]} # Record the numbers of reactions in each system
    sim_properties['rxns_numbers'] = rxns_numbers
    
    return None

#########################################################################################


#########################################################################################
def getlocusNumtoGeneSeq(sim_properties, gbfile):
    """
    
    Called at the very beginning of the simulation once

    Description: Set up sub dictionary sim_properties['locusNumtoIndex'] where the values are the genes' locusNums and keys are the start and end position of each gene
    """

    dna = gbfile
    gene_secondhalf = []
    gene_firsthalf = []
    locusNumtoIndex = {}
    locusNumtoGeneSeq = {}

    endposition = 543086
    # the index of the end position of 0910 is 543086
    # The length of the whole chromosome is len(dna.seq) is 543379
    position = 0

    gene_list = []
    for i in range(len(dna.features)):
        if ('product' in dna.features[i].qualifiers.keys()):
            #print(i) # This first statement works
            #print(dna.features[i].qualifiers['product'])
            if dna.features[i].qualifiers['product'][0]:# Figure out how to sort out for ribosomal operons?
                #print(dna.features[i].qualifiers['product'])
                gene_list.append(i)

    for gene in gene_list:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real
        if start < len(dna.seq)/2:
            gene_firsthalf.append(gene)
            if start == 0:
                locusNumtoIndex[locusNum] = [position, end]
                position = end
            
            else:
                locusNumtoIndex[locusNum] = [position, end]
                position = end
        else:
            gene_secondhalf.append(gene)

    position = endposition
    gene_secondhalf.reverse()

    for gene in gene_secondhalf:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real
        if end == endposition:

            locusNumtoIndex[locusNum] = [start, position]
            position = start
            
        else:
            locusNumtoIndex[locusNum] = [start, position]
            position = start

    sim_properties['locusNumtoIndex'] = locusNumtoIndex

    for locusNum, Index in locusNumtoIndex.items():
        start, end = Index[0], Index[1]
        gene_seq = Seq(str(dna.seq[start:end]))
        locusNumtoGeneSeq[locusNum] = gene_seq
    
    sim_properties['locusNumtoGeneSeq'] = locusNumtoGeneSeq

    return None

#########################################################################################




#########################################################################################
def gettRNAMap(sim_properties):
    """
    Input: sim_properties

    Return: None

    Called at the beginning of the simulation

    Description: sim_properties['trna_map']
    Set up the sub dictionary with amino acids and correspoding tRNAs {'LEU': ['R_0070', 'R_0423', 'R_0506'],...}
    """

    genome = sim_properties['genome']

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
                
    sim_properties['trna_map'] = tRNA_map

    return None

#########################################################################################

#########################################################################################
def getriboTolocusMap(sim_properties):
    """
    Input: 
    
    Description: mapping ribosomal protein names, e.g. L27 to their locusTags
    L33 has two genes, 0930 and 0932

    """
    
    ribotoLocus = {}

    for locusTag, subdic in sim_properties['genome'].items():
        if subdic['Type'] != 'gene':
            GeneName = subdic['GeneName']
            if GeneName.find('ribosomal protein') > 0:
                ribo = GeneName.split(' ')[-1]
                ribotoLocus[ribo] = []

    for locusTag, subdic in sim_properties['genome'].items():
        if subdic['Type'] != 'gene':
            GeneName = subdic['GeneName']
            if GeneName.find('ribosomal protein') > 0:
                ribo = GeneName.split(' ')[-1]
                locusNum = locusTag.split('_')[1]
                ribotoLocus[ribo].append(locusNum)
    
    # Manually Correct L27 since P_500 with name Maturation protease for ribosomal protein L27
    ribotoLocus['L27'] = ['0499']

    sim_properties['riboToLocus'] = ribotoLocus

    rPtn_locusNums = [locusNum for locusNums in ribotoLocus.values() for locusNum in locusNums]
    
    sim_properties['rPtn_locusNums'] = rPtn_locusNums


    return None


#########################################################################################

#########################################################################################
def getComplexMap(sim_properties):
    """
    Input: 
        cplx Excel Sheet

    Initialize the complex dictionary: {'ATPSynthase': {'Stoi':{'0789':1, ...},
    'init_count': 0, 'MP_count': 13:}, ... }
    
    """

    transmemPtnList = sim_properties['locations_ptns']['trans-membrane']
    cplx_path =  sim_properties['cplx_path']

    cplx_df = pd.read_excel(cplx_path, sheet_name='Complexes', dtype=str)
    # abc_df = pd.read_excel(cplx_path, sheet_name='ABC Transporter', dtype=str)
    # ATP_df = pd.read_excel(cplx_path, sheet_name='ATPSynthase', dtype=str)

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
        
        for locusNum in set(locusNums):
            cplx_dict[name]['Stoi'][locusNum] = cplx_composition.count(locusNum)


    for name, subdict in cplx_dict.items():
        Stoi = subdict['Stoi']
        transMP_count = 0
        for locusNum, stoi in Stoi.items():
            if locusNum in transmemPtnList:
                transMP_count += stoi
        subdict['transMP_count'] = transMP_count

    sim_properties['cplx_dict'] = cplx_dict
    
    print(f"cplx_dict {len(cplx_dict)}")

    for cplx, subdict in cplx_dict.items():
        print(cplx, subdict)

    return None
#########################################################################################

#########################################################################################
def outPtnCorrInitCounts(sim_properties, rank):
    """
    Description: Output the corrected counts for proteins in complexes
    """
    if rank == 1: # Only output the Excel sheet from one replicate

        cplx_dict = sim_properties['cplx_dict']

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
        
        def getPtnDict(row, locusNum, corrected_PtnIniCount, total_PtnIniCount, complex, AASequence):
            ptn = {}

            ptn['Locus Num'] = locusNum
            ptn['Gene Name'] = row['Gene Name']
            ptn['Gene Product'] = row['Gene Product']
            ptn['Localization'] = row['Localization']
            ptn['Protein Length'] = len(AASequence)

            ptn['Complex'] = ', '.join(a for a in complex)

            cplx_count = getCplxcounts(complex)
            ptn['Sim. Initial Cplx Cnt'] = ', '.join(str(a) for a in cplx_count)

            Stois = getStois(complex, locusNum)
            ptn['Stois in Cplx'] = ', '.join(str(a) for a in Stois)

            ptn['w/o. Cplx Sim. Ptn Cnt'] = row['Sim. Initial Ptn Cnt']

            ptn['w/. Cplx Sim. Free Ptn Cnt'] = corrected_PtnIniCount
            ptn['w. Cplx Sim. Total Ptn Cnt'] = total_PtnIniCount

            return ptn

        proteomicsDF = pd.read_excel(sim_properties['init_conc_path'], sheet_name='Comparative Proteomics', skiprows=[1])

        ptns = []

        for index, row in proteomicsDF.iterrows():
            locusTag = row['Locus Tag']; locusNum = locusTag.split('_')[1]

            PtnIniCount = row['Sim. Initial Ptn Cnt'] # number without complex considered
            corrected_PtnIniCount, total_PtnIniCount, complex = rxns_CME.correctInitPtnCount(PtnIniCount, locusNum, sim_properties)

            AASequence = sim_properties['genome'][locusTag]['AAsequence']
            ptn = getPtnDict(row, locusNum, corrected_PtnIniCount, total_PtnIniCount, complex, AASequence)

            ptns.append(ptn)


        sheet_names = ['Corrected Ptns Init Cnt']
        dfs = []
        dfs.append(pd.DataFrame(ptns).sort_values(by='Locus Num'))

        with pd.ExcelWriter('./Corrected_initial_ptns_count.xlsx') as writer:
            for sheet_name, df in zip(sheet_names, dfs):
                df.to_excel(writer, sheet_name=sheet_name, index=False)

    return None
#########################################################################################


#########################################################################################
def getTransloconRNCComplex(sim_properties):
    """
    Description: Get the dictionary of complexes with translocon Sec or RNaseY Dimer involved
    """

    Translocon_RNC_dict = {}

    transmemPtnList = sim_properties['locations_ptns']['trans-membrane']
    for locusNum in transmemPtnList:
        if locusNum not in sim_properties['sole_YidC_ptns']:
            SecYEGDF_RNC = 'SecYEGDF_RNC_'+locusNum
            Translocon_RNC_dict[SecYEGDF_RNC] = {}
            Translocon_RNC_dict[SecYEGDF_RNC]['init_count'] = 0
            Translocon_RNC_dict[SecYEGDF_RNC]['Stoi'] = {'P_0652':1, 'P_0839':1, 'P_0774':1, 'P_0412':1, 'Ribosome':1}
            Translocon_RNC_dict[SecYEGDF_RNC]['transMP_count'] = 4

            Sec_SR_SRP_RNC = 'SecYEGDF_SR_SRP_RNC_' + locusNum
            Translocon_RNC_dict[Sec_SR_SRP_RNC] = {}
            Translocon_RNC_dict[Sec_SR_SRP_RNC]['init_count'] = 0
            Translocon_RNC_dict[Sec_SR_SRP_RNC]['Stoi'] = {'P_0652':1, 'P_0839':1, 'P_0774':1, 'P_0412':1, 'Ribosome':1, 'P_0360':1, 'P_0429':1}
            Translocon_RNC_dict[Sec_SR_SRP_RNC]['transMP_count'] = 4

            
        else:
            # YidC_RNC_id = 'YidC_RNC_' + locusNum
            YidC_CPtn = 'YidC_CP_'+locusNum
            Translocon_RNC_dict[YidC_CPtn] = {}
            Translocon_RNC_dict[YidC_CPtn]['init_count'] = 0
            Translocon_RNC_dict[YidC_CPtn]['Stoi'] = {'P_0908':1}
            Translocon_RNC_dict[YidC_CPtn]['transMP_count'] = 1

    lipoproteinList = sim_properties['locations_ptns']['lipoprotein']
    secretedList = sim_properties['locations_ptns']['extracellular']

    for locusNum in lipoproteinList+secretedList:

        SecYEGDF_SecA_RNC = 'SecYEGDF_SecA_RNC_' + locusNum
        Translocon_RNC_dict[SecYEGDF_SecA_RNC] = {}
        Translocon_RNC_dict[SecYEGDF_SecA_RNC]['init_count'] = 0
        Translocon_RNC_dict[SecYEGDF_SecA_RNC]['Stoi'] = {'P_0652':1, 'P_0839':1, 'P_0774':1, 'P_0412':1,'Ribosome':1, 'P_0095':1}
        Translocon_RNC_dict[SecYEGDF_SecA_RNC]['transMP_count'] = 4

    for locusNum in lipoproteinList:
        Lgt1_PreP = 'Lgt1_PreP_' + locusNum
        Lgt2_PreP = 'Lgt2_PreP_' + locusNum
        Lsp_ProP = 'Lsp_ProP_' + locusNum 
        Translocon_RNC_dict[Lgt1_PreP] = {}
        Translocon_RNC_dict[Lgt1_PreP]['init_count'] = 0
        Translocon_RNC_dict[Lgt1_PreP]['Stoi'] = {'P_0818':1, f'lipo_{locusNum}':1}
        Translocon_RNC_dict[Lgt1_PreP]['transMP_count'] = 1

        Translocon_RNC_dict[Lgt2_PreP] = {}
        Translocon_RNC_dict[Lgt2_PreP]['init_count'] = 0
        Translocon_RNC_dict[Lgt2_PreP]['Stoi'] = {'P_0820':1, f'lipo_{locusNum}':1}
        Translocon_RNC_dict[Lgt2_PreP]['transMP_count'] = 1

        Translocon_RNC_dict[Lsp_ProP] = {}
        Translocon_RNC_dict[Lsp_ProP]['init_count'] = 0
        Translocon_RNC_dict[Lsp_ProP]['Stoi'] = {'P_0518':1, f'lipo_{locusNum}':1}
        Translocon_RNC_dict[Lsp_ProP]['transMP_count'] = 1

        # SecYEGDF_SecA_CP = 'SecYEGDF_SecA_CP_' + locusNum
        # Translocon_RNC_dict[SecYEGDF_SecA_CP] = {}
        # Translocon_RNC_dict[SecYEGDF_SecA_CP]['init_count'] = 0
        # Translocon_RNC_dict[SecYEGDF_SecA_CP]['Stoi'] = {'P_0652':1, 'P_0839':1, 'P_0774':1, 'P_0412':1, 'P_0095':1}
        # Translocon_RNC_dict[SecYEGDF_SecA_CP]['transMP_count'] = 4

    for locusTag, locusDict in sim_properties['genome'].items():
        locusNum = locusTag.split('_')[1]
        if locusDict['Type'] == 'protein': # mRNA
            Degradosome_mRNA = 'Degradosome_mRNA_' + locusNum
            Translocon_RNC_dict[Degradosome_mRNA] = {}
            Translocon_RNC_dict[Degradosome_mRNA]['init_count'] = 0
            Translocon_RNC_dict[Degradosome_mRNA]['Stoi'] = {'P_0359':1, 'P_0600':1, 'P_257':1, 'P_437':1, f"mRNA_{locusNum}":1}
            Translocon_RNC_dict[Degradosome_mRNA]['transMP_count'] = 1


    # Proteins subject to FtsH protease

    FtsH_CP_list = set(transmemPtnList + lipoproteinList) - set(sim_properties['sole_YidC_ptns'])
    # print('FtsH',len(FtsH_CP_list), FtsH_CP_list)

    for locusNum in FtsH_CP_list:
        FtsH_CP = 'FtsH_CP_' + locusNum
        Translocon_RNC_dict[FtsH_CP] = {}
        Translocon_RNC_dict[FtsH_CP]['init_count'] = 0
        Translocon_RNC_dict[FtsH_CP]['Stoi'] = {'P_0039':1, f'P_{locusNum}':1}
        Translocon_RNC_dict[FtsH_CP]['transMP_count'] = 1



    sim_properties['Translocon_RNC_dict'] = Translocon_RNC_dict
    

    return None
#########################################################################################


#########################################################################################
def initializeCosts(sim_properties):
    """
    
    Called once at the beginning of the simulation
    
    Description: Initialize the nucleotides cost species in counts dictionary

    """

    countsDic = sim_properties['counts']

    # nucleotides cost as energetic reactions NTPs -> NDPs + pi; For tRNA charging to AMP and ppi
    energy_cost_counters = {'GTP_translate_cost':'g', 'GTP_transloc_cost':'g', 'ATP_trsc_cost':'a', 'ATP_mRNAdeg_cost':'a', 
                            'ATP_DNArep_cost':'a', 'ATP_transloc_cost':'a','ATP_tRNAcharging_cost':'a', 'ATP_Ptndeg_cost':'a'}

    # nucleotides cost in polymerization processes NTPs -> NMPs + ppi
    nuc_costs = {'ATP_mRNA_cost':'M_atp_c', 'CTP_mRNA_cost':'M_ctp_c', 'UTP_mRNA_cost':'M_utp_c', 'GTP_mRNA_cost':'M_gtp_c',
                    'ATP_tRNA_cost':'M_atp_c', 'CTP_tRNA_cost':'M_ctp_c', 'UTP_tRNA_cost':'M_utp_c', 'GTP_tRNA_cost':'M_gtp_c',
                    'ATP_rRNA_cost':'M_atp_c', 'CTP_rRNA_cost':'M_ctp_c', 'UTP_rRNA_cost':'M_utp_c', 'GTP_rRNA_cost':'M_gtp_c',
                    'dATP_DNArep_cost':'M_datp_c', 'dTTP_DNArep_cost':'M_dttp_c', 'dCTP_DNArep_cost':'M_dctp_c', 'dGTP_DNArep_cost':'M_dgtp_c'}
    #
    NMP_recycle_counters = {'AMP_mRNAdeg_recycled':'M_amp_c', 'UMP_mRNAdeg_recycled':'M_ump_c', 'CMP_mRNAdeg_recycled':'M_cmp_c', 'GMP_mRNAdeg_recycled':'M_gmp_c'}


    for key, value in energy_cost_counters.items():
        countsDic[key] = [int(0)]
        # countsDic[value+'_accumulative'] = [int(0)]
    
    for key, value in nuc_costs.items():
        countsDic[key] = [int(0)]
        # countsDic[value+'_accumulative'] = [int(0)]

    for key, value in NMP_recycle_counters.items():
        countsDic[key] = [int(0)]
        # countsDic[value+'_accumulative'] = [int(0)]
    
    monomerlist = ['M_atp_c', 'M_utp_c', 'M_ctp_c', 'M_gtp_c','M_datp_c','M_dttp_c','M_dctp_c','M_dgtp_c']
    
    for metID in monomerlist:
        # nucName = monomerName.split('_')[1]
        NTP_shortage = metID + '_shortage'
        countsDic[NTP_shortage] = [int(0)]
    
    # Amino acids
    letter2aa = {'A': 'ALA','R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 
                 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 
                 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
    aaCostMap = {}; unpaidaaCostMap = {}; recycledaaMap = {}; lettertoAAMet = {}

    for letter, aa in  letter2aa.items():
        aaCostMap[letter] = aa + '_cost'
        unpaidaaCostMap[letter] = aa + '_cost_unpaid'
        recycledaaMap[letter] = aa + '_ptndeg_recycled'
        if aa == 'GLY':
            lettertoAAMet[letter] = 'M_***_c'.replace('***', aa.lower())
        else:
            lettertoAAMet[letter] = 'M_***__L_c'.replace('***', aa.lower())

    sim_properties['letter2aa'] = letter2aa
    sim_properties['aaCostMap'] = aaCostMap
    sim_properties['unpaidaaCostMap'] = unpaidaaCostMap
    sim_properties['recycledaaMap'] = recycledaaMap
    sim_properties['lettertoAAMet'] = lettertoAAMet

    for aa_recycled in sim_properties['recycledaaMap'].values():
        countsDic[aa_recycled] = [int(0)]

    # sim_properties['aaCostMap'] = {"A":"ALA_cost", "R":"ARG_cost", 
    #     "N":"ASN_cost", "D":"ASP_cost", "C":"CYS_cost", "E":"GLU_cost", "Q":"GLN_cost", "G":"GLY_cost", 
    #         "H":"HIS_cost", "I":"ILE_cost", "L":"LEU_cost", "K":"LYS_cost", "M":"MET_cost", "F":"PHE_cost", 
    #     "P":"PRO_cost", "S":"SER_cost", "T":"THR_cost", "W":"TRP_cost", "Y":"TYR_cost", "V":"VAL_cost"}

    # # Mapping aa sequence to unpaid aacost species in CME
    # sim_properties['unpaidaaCostMap'] = {'A': 'ALA_cost_unpaid', 'R': 'ARG_cost_unpaid', 'N': 'ASN_cost_unpaid', 'D': 'ASP_cost_unpaid', 'C': 'CYS_cost_unpaid',
    #   'E': 'GLU_cost_unpaid', 'Q': 'GLN_cost_unpaid', 'G': 'GLY_cost_unpaid', 'H': 'HIS_cost_unpaid', 'I': 'ILE_cost_unpaid', 
    #   'L': 'LEU_cost_unpaid', 'K': 'LYS_cost_unpaid', 'M': 'MET_cost_unpaid', 'F': 'PHE_cost_unpaid', 'P': 'PRO_cost_unpaid', 
    #   'S': 'SER_cost_unpaid','T': 'THR_cost_unpaid', 'W': 'TRP_cost_unpaid', 'Y': 'TYR_cost_unpaid', 'V': 'VAL_cost_unpaid'}



    return None 
#########################################################################################



#########################################################################################
def initializeMediumConcs(sim_properties):
    
    sim_properties['medium'] = {}
    
    sim_medium = pd.read_excel(sim_properties['init_conc_path'], sheet_name='Simulation Medium')
    
    for row, nutrient in sim_medium.iterrows():
        
        metID = 'M_' + nutrient['Met ID']
        
        # The concentration of medium species is constant and put in another subDic.
        sim_properties['medium'][metID] = nutrient['Conc (mM)']
        

    return None
#########################################################################################


#########################################################################################
def initializeMetabolitesCounts(sim_properties):
    """
    Input: sim_properties

    Return: None

    Description: Add the initial counts of metabolites into counts dictionary
    
    """


    metabolite_ic = pd.read_excel(sim_properties['init_conc_path'], sheet_name='Intracellular Metabolites')
    countsDic = sim_properties['counts']

    for row, metabolite in metabolite_ic.iterrows():
        
        metID = 'M_' + metabolite['Met ID']
        
        # The first element of metabolite's trajectory is the initial counts at 0 second
        countsDic[metID] = [mMtoPart(metabolite['Init Conc (mM)'], sim_properties)]
    
    # print('Metabolism Initialized')
    return None
#########################################################################################


#########################################################################################
def initializeProteinMetabolitesCounts(sim_properties):
    """

    Called once at the very beginning of the simulation

    Description: Define and add the initial counts of forms of four proteins 'P_0233', 'P_0694', 'P_0234', 'P_0779' in phosphorelay; The phospholated form starts from 0
    """

    phosphorelayPtns = ['P_0233', 'P_0694', 'P_0234', 'P_0779']

    phosphorelayPtnsMap = {'P_0233':352, 'P_0694':289, 'P_0234':313, 'P_0779':830}

    # phosphorelayPtnsCounts = [352, 289, 313, 830]

    data_file =  sim_properties['init_conc_path']
    
    ptnMets = pd.read_excel(data_file, sheet_name='Protein Metabolites')
    
    for index, row in ptnMets.iterrows():
        PtnID = row['Protein']
        if PtnID in phosphorelayPtns:

            metabolites = row['Metabolite IDs'].split(',')
        
            for metID in metabolites:
                if metID == metabolites[0]:

                    sim_properties['counts'][metID] = [phosphorelayPtnsMap[PtnID]]

                else:

                    sim_properties['counts'][metID] = [0]
        

    return None
#########################################################################################

#########################################################################################
def getGIPSpeicesCounts(sim, sim_properties):
    """

    Description: Combine and remove the duplicates in the lists of species and counts when adding GIP reactions 
                Compare the recored species list with sim.particleMap
    """

    pre_list = sim_properties['pre_list']; pre_counts = sim_properties['pre_counts']
    ini_list = sim_properties['ini_list']; ini_counts = sim_properties['ini_counts'] 
    rep_list = sim_properties['rep_list']; rep_counts = sim_properties['rep_counts'] 
    trsc_list = sim_properties['trsc_list']; trsc_counts = sim_properties['trsc_counts']
    translation_list = sim_properties['translation_list']; translation_counts = sim_properties['translation_counts']
    Deg_list = sim_properties['Deg_list']; Deg_counts = sim_properties['Deg_counts']
    tRNA_list = sim_properties['tRNA_list']; tRNA_counts = sim_properties['tRNA_counts']
    ribo_list = sim_properties['ribo_list']; ribo_counts = sim_properties['ribo_counts']
    cplx_list = sim_properties['cplx_list']; cplx_counts = sim_properties['cplx_counts']

    # lists = [pre_list, ini_list, rep_list, trsc_list, translation_list, Deg_list, tRNA_list, ribo_list, cplx_list]
    # counts = [pre_counts, ini_counts, rep_counts, trsc_counts, translation_counts, Deg_counts, tRNA_counts, ribo_counts, cplx_counts]

    entire_species = list(chain(pre_list, ini_list, rep_list, trsc_list, translation_list, Deg_list, tRNA_list, ribo_list, cplx_list))
    entire_counts = list(chain(pre_counts, ini_counts, rep_counts, trsc_counts, translation_counts, Deg_counts, tRNA_counts, ribo_counts, cplx_counts))

    GIP_species = []; GIP_initial_counts = []
    seen = set()
    for specie, count in zip(entire_species, entire_counts):
        if specie not in seen:
            seen.add(specie)
            GIP_species.append(specie)
            GIP_initial_counts.append(count)

    defined_species = list(sim.particleMap.keys())

    if set(GIP_species) - set(defined_species) != set():
        not_defined_species = list(set(GIP_species)- set(defined_species))
        print(f"Error: {len(not_defined_species)} were not defined in CME simulation: {','.join(not_defined_species)}")
    elif set(defined_species) - set(GIP_species) != set():
        not_recored_species = list(set(defined_species)- set(GIP_species))
        print(f"Error: {len(not_recored_species)} were not recored into sim_properties: {','.join(not_recored_species)}")
    else:
        print(f'CHECKPOINT: Total {len(GIP_species)} species defined in CME simulation and added into sim_properties') 

    # duplicated_species = list(set(entire_species) - set(GIP_species))

    sim_properties['GIP_species'] = GIP_species
    sim_properties['GIP_initial_counts'] = GIP_initial_counts

    return GIP_species, GIP_initial_counts
#########################################################################################


#########################################################################################
def addGeneticInformationSpeciesCounts(sim,sim_properties):
    """
    Input: sim, sim_properties

    Return: None

    Called when restart new CME simulation

    Description: Add species counts to new CME simulation; Initialize the trajectory of species in sim_properties dictionary; 
    """


    GIP_species, GIP_initial_counts = getGIPSpeicesCounts(sim, sim_properties)

    time_second = sim_properties['time_second']
    countsDic = sim_properties['counts']

    if time_second[-1] == 0:
        for i in range(len(GIP_species)):
            sim.addParticles(species = GIP_species[i], count = GIP_initial_counts[i])
            
            countsDic[GIP_species[i]] = [GIP_initial_counts[i]]

    else:
        for i in range(len(GIP_species)):
            sim.addParticles(species = GIP_species[i], count = countsDic[GIP_species[i]][-1])
    
    # addpreCounts(sim, sim_properties)
    # # addEnzymesCounts(sim, sim_properties)
    # addinitiationCounts(sim, sim_properties)
    # addReplicationCounts(sim, sim_properties)
    # addTranscriptionCounts(sim, sim_properties)
    # addTranslationCounts(sim, sim_properties)
    # addDegradationCounts(sim, sim_properties)
    # addtRNAChargingCounts(sim, sim_properties)
    # # addaaCostCounts(sim, sim_properties)
    # addriboCounts(sim, sim_properties)
    # addcplxCounts(sim, sim_properties)
    
    return None

#########################################################################################


# #########################################################################################
# def addinitiationCounts(sim, sim_properties):
#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     ini_list = sim_properties['ini_list']
#     ini_counts = sim_properties['ini_counts']

#     if time_second[-1] == 0:
#         for i in range(len(ini_list)):
#             sim.addParticles(species = ini_list[i], count = ini_counts[i])
#             countsDic[ini_list[i]] = [ini_counts[i]]

#     else:
#         for i in range(len(ini_list)):
#             sim.addParticles(species = ini_list[i], count = countsDic[ini_list[i]][-1])

#     return None
# #########################################################################################



# #########################################################################################
# def addReplicationCounts(sim, sim_properties):
#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     rep_list = sim_properties['rep_list']
#     rep_counts = sim_properties['rep_counts']


#     if time_second[-1] == 0:
#         for i in range(len(rep_list)):
#             sim.addParticles(species = rep_list[i], count = rep_counts[i])
#             countsDic[rep_list[i]] = [rep_counts[i]]

#     else:
#         for i in range(len(rep_list)):
#             sim.addParticles(species = rep_list[i], count = countsDic[rep_list[i]][-1])
#     return None
# #########################################################################################


# #########################################################################################
# def addTranscriptionCounts(sim, sim_properties):

#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']

#     trsc_list = sim_properties['trsc_list']
#     trsc_counts = sim_properties['trsc_counts']


#     if time_second[-1] == 0:

#         for i in range(len(trsc_list)):
#             # print(trsc_list[i])
#             sim.addParticles(species = trsc_list[i], count = trsc_counts[i])
#             countsDic[trsc_list[i]] = [trsc_counts[i]]
        
#         # print("Transcription Initialized")

#     else:
#         for i in range(len(trsc_list)):
#             sim.addParticles(species = trsc_list[i], count = countsDic[trsc_list[i]][-1])

#     return None
# #########################################################################################



# #########################################################################################
# def addTranslationCounts(sim, sim_properties):
    
#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']

#     translation_list = sim_properties['translation_list']
#     translation_counts = sim_properties['translation_counts']
#     # When defining the number of each species, we need to give the a conditional statement 
#     # to differ the very beginning of the whole cycle and the restart simulation
#     if time_second[-1] == 0:

#         for i in range(len(translation_list)):
#             sim.addParticles(species = translation_list[i], count = translation_counts[i])
#             countsDic[translation_list[i]] = [translation_counts[i]]
#         # print('Translation Initialized')
#     else:
        
#         # Pass the last element of trajectory recorded in coundsDic to next newly restart CME simulation
#         for i in range(len(translation_list)):
#             sim.addParticles(species = translation_list[i], count = countsDic[translation_list[i]][-1])

#     #print('Translation Initialized')

#     return None

# #########################################################################################



# #########################################################################################
# def addDegradationCounts(sim, sim_properties):
#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     Deg_list = sim_properties['Deg_list']
#     Deg_counts = sim_properties['Deg_counts']

#     if time_second[-1] == 0:
#         for i in range(len(Deg_list)):
#             sim.addParticles(species = Deg_list[i], count = Deg_counts[i])
#             countsDic[Deg_list[i]] = [Deg_counts[i]]
#     else:
#         for i in range(len(Deg_list)):
#             sim.addParticles(species = Deg_list[i], count = countsDic[Deg_list[i]][-1])


#     return None

# #########################################################################################

# #########################################################################################
# def addtRNAChargingCounts(sim, sim_properties):
#     """
    
    
#     Description: addParticles to tRNA charging related species, including atp, amp, ppi, R_XXXX_ch
#     """
    
#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     tRNA_list = sim_properties['tRNA_list'] 
#     tRNA_counts = sim_properties['tRNA_counts']
    
#     if time_second[-1] == 0:
#         for i in range(len(tRNA_list)):
#             sim.addParticles(species = tRNA_list[i], count = tRNA_counts[i])
#             countsDic[tRNA_list[i]] = [tRNA_counts[i]]
#     else:
#         for i in range(len(tRNA_list)):
#             sim.addParticles(species = tRNA_list[i], count = countsDic[tRNA_list[i]][-1])
#     return None

# #########################################################################################

# def addaaCostCounts(sim,sim_properties):
#     """
    
#     Description: addParticles to aa_cost, aa_cost_unpaid species to CME; initialize aa_cost, aa_cost_unpaid in countsDic
#     """

#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     aaCost_list = sim_properties['aaCost_list']
#     unpaidaaCost_list = sim_properties['unpaidaaCost_list']

#     if time_second[-1] == 0:
#         for i in range(len(aaCost_list)):
#             sim.addParticles(species = aaCost_list[i], count = 0)
#             countsDic[aaCost_list[i]] = [0]
#         for i in range(len(unpaidaaCost_list)):
#             sim.addParticles(species = unpaidaaCost_list[i], count = 0)
#             countsDic[unpaidaaCost_list[i]] = [0]

#     else:
#         for i in range(len(aaCost_list)):
#             sim.addParticles(species = aaCost_list[i], count = countsDic[aaCost_list[i]][-1])
   
#         for i in range(len(unpaidaaCost_list)):
#             sim.addParticles(species = unpaidaaCost_list[i], count = countsDic[unpaidaaCost_list[i]][-1])

#     return None


# #########################################################################################
# def addriboCounts(sim, sim_properties):
#     """
#     Description: add intermediates in ribosome biogenesis
    
#     """

#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     ribo_list = sim_properties['ribo_list']
#     ribo_counts = sim_properties['ribo_counts']
#     if time_second[-1] == 0:
#         for i in range(len(ribo_list)):
#             sim.addParticles(species = ribo_list[i], count = ribo_counts[i])
#             countsDic[ribo_list[i]] = [ribo_counts[i]]
#     else:
#         for i in range(len(ribo_list)):
#             sim.addParticles(species = ribo_list[i], count = countsDic[ribo_list[i]][-1])

#     # print('Counts of species in ribosome biogenesis added at time {0}'.format(time_second[-1]))

#     return None


# #########################################################################################

# #########################################################################################
# def addcplxCounts(sim, sim_properties):

#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     cplx_list = sim_properties['cplx_list']
#     cplx_counts = sim_properties['cplx_counts']

#     if time_second[-1] == 0:
#         for i in range(len(cplx_list)):
#             sim.addParticles(species = cplx_list[i], count = cplx_counts[i])
#             countsDic[cplx_list[i]] = [cplx_counts[i]]
#     else:
#         for i in range(len(cplx_list)):
#             sim.addParticles(species = cplx_list[i], count = countsDic[cplx_list[i]][-1])

#     print('Counts of species in protein complexes added at time {0}'.format(time_second[-1]))

#     return None

# #########################################################################################

# #########################################################################################
# def addpreCounts(sim, sim_properties):
#     """
#     Description: add counts to predefined species
    
#     """

#     time_second = sim_properties['time_second']
#     countsDic = sim_properties['counts']
#     pre_list = sim_properties['pre_list']
#     pre_counts = sim_properties['pre_counts']
#     if time_second[-1] == 0:
#         for i in range(len(pre_list)):
#             sim.addParticles(species = pre_list[i], count = pre_counts[i])
#             countsDic[pre_list[i]] = [pre_counts[i]]
#     else:
#         for i in range(len(pre_list)):
#             sim.addParticles(species = pre_list[i], count = countsDic[pre_list[i]][-1])

#     # print('Counts of species in ribosome biogenesis added at time {0}'.format(time_second[-1]))

#     return None
# #########################################################################################


# #########################################################################################
# def addEnzymesCounts(sim, sim_properties):
#     # Initialize the numbers of three enzymes when new CME simulation is created
#     countsDic = sim_properties['counts']

    
#     if sim_properties['time_second'][-1] == 0:
#         occupied_complexes = sim_properties['occupied_complexes']
#         print('Occupied_complexes: {0}'.format(occupied_complexes))

#         occupied_RNAP = occupied_complexes['occupied_RNAP']
#         occupied_ribosome = occupied_complexes['occupied_ribosome']
#         occupied_degra = occupied_complexes['occupied_degradosome']

#         # Need to define species
#         # The numbers for the very beginning of whole cell cycle
#         sim.addParticles(species = 'RNAP', count = 187 - occupied_RNAP)
#         sim.addParticles(species = 'Ribosome', count = 500 - occupied_ribosome)
#         sim.addParticles(species = 'Degradosome', count = 120 - occupied_degra)

#         # For each species, a list is created with initial value
#         countsDic['RNAP'] = [int(187 - occupied_RNAP)]
#         countsDic['Ribosome'] = [int(500 - occupied_ribosome)]
#         countsDic['Degradosome'] = [int(120 - occupied_degra)]

#     else:
#         # Read the last element of the list and pass to next restart simulation
#         sim.addParticles(species = 'RNAP',count = countsDic['RNAP'][-1] )
#         sim.addParticles(species = 'Ribosome',count = countsDic['Ribosome'][-1] )
#         sim.addParticles(species = 'Degradosome',count = countsDic['Degradosome'][-1] )


# #########################################################################################




#########################################################################################
def mMtoPart(conc, sim_properties):
    """
    Convert particle counts to mM concentrations for the ODE Solver

    Parameters:
    particles (int): The number of particles for a given chemical species

    Returns:
    conc (float): The concentration of the chemical species in mM
    """

    ### Constants
    NA = 6.022e23 # Avogadro's
    
    count = max(1, int(round((conc/1000)*NA*sim_properties['volume_L'][-1])))

    return count
#########################################################################################




#########################################################################################
def initializeMembrane(sim_properties):
    """
    Input: sim_properties dictionary

    Return: None

    Description: 
    initialize the sim_properties['SA'] subdictionary to record the SA and volume traces
    Calculate the lipids' SA contribution and average membrane protein subunit SA contribution
    Calcualte the cellVolume based on 200 nm radius
    
    """
    countsDic = sim_properties['counts']

    # avgProtSA = 28 #24.75 # nm^2, average protein surface area to produce expected 47% coverage for 9.6K membrane proteins
    
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
    
    lipid_area = 0
    
    for lipid, size in lipidSizes.items():
        
        lipid_area = lipid_area + countsDic[lipid][-1]*size
    
    # Double layered lipids
    lipidSA = int(lipid_area*0.513) #nm^2

    # Calculate the surface area at the initial time
    if sim_properties['time_second'][-1] == 0:
        avg_MP_area = getAvgPtnArea(sim_properties, lipidSA)

    avg_MP_area = sim_properties['avg_MP_area']

    MP_total_count, MP_cplx_count, MP_subcplx_count, MP_translocon_RNC, MP_free_count = getMPCounts(sim_properties)

    MPSA = int(avg_MP_area*MP_total_count) # nm^2

    sim_properties['SA'] = {}

    sim_properties['SA']['SA_lipid'] = [lipidSA]

    sim_properties['SA']['SA_ptn'] = [MPSA]
    
    sim_properties['SA']['SA_nm2'] = [int(lipidSA + MPSA)]
    
    sim_properties['SA']['SA_m2'] = [int(lipidSA + MPSA )/1e18]    # m^2

    sim_properties['division_started'] = [False]

    return None


def getAvgPtnArea(sim_properties, lipid_SA):
    """
    Description: Calculate the average surface area contribution per membrane protein/subunit at the beginning
    """
    
    MP_total_count, MP_cplx_count, MP_subcplx_count, MP_translocon_RNC, MP_free_count = getMPCounts(sim_properties)

    print(f"{MP_cplx_count} transmembrane protein subunits assembly into final complex, {MP_subcplx_count} into sub complex, {MP_translocon_RNC} into Translocon_RNC complex, {MP_free_count} are free")
    
    MP_total_count = MP_cplx_count + MP_subcplx_count + MP_translocon_RNC + MP_free_count

    print(f"Total {MP_total_count} transmembrane protein at begenning")

    r_cell_nm = sim_properties['r_cell']*10**9

    SA = 4*np.pi*(r_cell_nm)**2
    MP_SA = SA - lipid_SA #nm^2

    print(f"Lipids account for {100*lipid_SA/SA:.2f} % percentage of Initial Surface Area")
    print(f"Transmembrane Proteins account for {100*MP_SA/SA:.2f} % percentage of Initial Surface Area")

    avg_MP_area = MP_SA/MP_total_count

    sim_properties['avg_MP_area'] = avg_MP_area

    print(f"Average surface area contribution per transmembrane protein is {avg_MP_area:.2f} nm^2 ")
    
    return avg_MP_area


def getMPCounts(sim_properties):
    """
    MP counts in final complex, subcomplex, and free
    """

    countsDic = sim_properties['counts']

    cplx_dict = sim_properties['cplx_dict']
    sub_cplx_dict = sim_properties['sub_cplx_dict']
    Translocon_RNC_dict = sim_properties['Translocon_RNC_dict']

    transmemPtnList = sim_properties['locations_ptns']['trans-membrane']

    MP_cplx_count = 0; MP_subcplx_count = 0; MP_translocon_RNC = 0; MP_free_count = 0

    for cplx, subdict in cplx_dict.items():
        transMP_count = subdict['transMP_count']
        # if sim_properties['time_second'][-1] != 0:
        #     if countsDic[cplx][-2] != countsDic[cplx][-1] and transMP_count != 0 :
        #         print(f"Time {sim_properties['time_second'][-1]}")
        #         print(f"{cplx} {countsDic[cplx][-2]} {countsDic[cplx][-1]}")

        MP_cplx_count += transMP_count*countsDic[cplx][-1]

    for subcplx, subdict in sub_cplx_dict.items():
        transMP_count = subdict['transMP_count']
        # if sim_properties['time_second'][-1] != 0:
        #     if countsDic[subcplx][-2] != countsDic[subcplx][-1] and transMP_count != 0 :
        #         print(f"Time {sim_properties['time_second'][-1]}")
        #         print(f"{subcplx} {countsDic[subcplx][-2]} {countsDic[subcplx][-1]}")

        MP_subcplx_count += transMP_count*countsDic[subcplx][-1]   
    
    for translocon_RNC, subdict in Translocon_RNC_dict.items():
        transMP_count = subdict['transMP_count'] 
        # if sim_properties['time_second'][-1] != 0:
        #     if countsDic[translocon_RNC][-2] != countsDic[translocon_RNC][-1] and transMP_count != 0 :
        #         print(f"Time {sim_properties['time_second'][-1]}")
        #         print(f"{translocon_RNC} {countsDic[translocon_RNC][-2]} {countsDic[translocon_RNC][-1]}")

        MP_translocon_RNC += transMP_count*countsDic[translocon_RNC][-1]   

    for MP_locusNum in transmemPtnList:
        ptnID = 'P_' + MP_locusNum
        MP_free_count += countsDic[ptnID][-1]

        # if sim_properties['time_second'][-1] != 0:
        #     if countsDic[ptnID][-2] != countsDic[ptnID][-1]:
        #         print(f"Time {sim_properties['time_second'][-1]}")
        #         print(f"{ptnID} {countsDic[ptnID][-2]} {countsDic[ptnID][-1]}")

    MP_total_count = MP_cplx_count + MP_subcplx_count + MP_translocon_RNC + MP_free_count

    print(f"At time {sim_properties['time_second'][-1]} second total MP count {MP_total_count}")
    print(f"{MP_cplx_count} transmembrane protein subunits assembly into final complex, {MP_subcplx_count} into sub complex, {MP_translocon_RNC} into Translocon_RNC complex, {MP_free_count} are free")

    # print(f"At time {sim_properties['time_second'][-1]} cplx MP count {MP_cplx_count}")
    # print(f"At time {sim_properties['time_second'][-1]} subcplx MP count {MP_subcplx_count}")
    # print(f"At time {sim_properties['time_second'][-1]} translocon RNC MP count {MP_translocon_RNC}")
    # print(f"At time {sim_properties['time_second'][-1]} free MP count {MP_free_count}")

    return MP_total_count, MP_cplx_count, MP_subcplx_count, MP_translocon_RNC, MP_free_count


# #########################################################################################
# def outPtnCorrInitCounts(sim_properties, rank):
#     """
#     Description: Output the corrected counts for proteins in complexes; Protein subunits in complex or not separated in two sheets
#     """
#     if rank == 1: # Only output the Excel sheet from one replicate

#         cplx_dict = sim_properties['cplx_dict']

#         def getCplxcounts(complex):
#             cplx_count = []

#             for cplx in complex:
#                 count = cplx_dict[cplx]['init_count']
#                 cplx_count.append(count)

#             return cplx_count

#         def getStois(complex, locusNum):
#             Stois = []

#             for cplx in complex:
#                 Stoi = cplx_dict[cplx]['Stoi'][locusNum]
#                 Stois.append(Stoi)

#             return Stois
        
#         def getPtnDict(row, locusNum, corrected_PtnIniCount, total_PtnIniCount, complex):
#             ptn = {}

#             if complex == []:
#                 ptn['Locus Num'] = locusNum
#                 ptn['Gene Name'] = row['Gene Name']
#                 ptn['Gene Product'] = row['Gene Product']
#                 ptn['Localization'] = row['Localization']
#                 ptn['Exp. Ptn Cnt'] = row['Exp. Ptn Cnt']
#                 ptn['Sim. Initial Ptn Cnt'] = corrected_PtnIniCount

#             else:
#                 ptn['Locus Num'] = locusNum
#                 ptn['Gene Name'] = row['Gene Name']
#                 ptn['Gene Product'] = row['Gene Product']
#                 ptn['Localization'] = row['Localization']

#                 ptn['Complex'] = ', '.join(a for a in complex)

#                 cplx_count = getCplxcounts(complex)
#                 ptn['Sim. Initial Cplx Cnt'] = ', '.join(str(a) for a in cplx_count)

#                 Stois = getStois(complex, locusNum)
#                 ptn['Stois in Cplx'] = ', '.join(str(a) for a in Stois)

#                 ptn['Exp. Ptn Cnt'] = row['Exp. Ptn Cnt']
#                 ptn['Sim. Free Initial Ptn Cnt'] = corrected_PtnIniCount
#                 ptn['Sim. Total Initial Ptn Cnt'] = total_PtnIniCount

#             return ptn

#         proteomicsDF = pd.read_excel(sim_properties['init_conc_path'], sheet_name='Comparative Proteomics', skiprows=[1])

#         changed_ptns = []
#         unchanged_ptns = []

#         for index, row in proteomicsDF.iterrows():
#             locusTag = row['Locus Tag']; locusNum = locusTag.split('_')[1]

#             PtnIniCount = row['Sim. Initial Ptn Cnt']
#             corrected_PtnIniCount, total_PtnIniCount, complex = rxns_CME.correctInitPtnCount(PtnIniCount, locusNum, sim_properties)

#             if complex == []:
#                 ptn = getPtnDict(row, locusNum, corrected_PtnIniCount, total_PtnIniCount, complex)
#                 unchanged_ptns.append(ptn)

#             else:
#                 ptn = getPtnDict(row, locusNum, corrected_PtnIniCount, total_PtnIniCount, complex)
#                 changed_ptns.append(ptn)

#         sheet_names = ['Corr Ptns Init Cnt', 'Ptns Init Cnt']
#         dfs = []
#         dfs.append(pd.DataFrame(changed_ptns).sort_values(by='Locus Num'))
#         dfs.append(pd.DataFrame(unchanged_ptns).sort_values(by='Locus Num'))

#         with pd.ExcelWriter('./Corrected_initial_ptns_count.xlsx') as writer:
#             for sheet_name, df in zip(sheet_names, dfs):
#                 df.to_excel(writer, sheet_name=sheet_name, index=False)

#     return None
# #########################################################################################