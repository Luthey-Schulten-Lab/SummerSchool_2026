"""
Source Code: rep_restart.py from CME_ODE model in Zane Thornburg, Cell 2022

Modification by Enguang Fu

Called by Rnxs_CME.py to add replication initiation and replication reactions when CME simulation initialized

Compared to the source code, 

    the unnecessary ssDNAboundSite_i species are removed.

    Update the rate constants based on instant volume

    Initiation finishes with dnaA filament length above 20

Zhang 2017 Paper:
    The presence of Mg2+ and ddATP can affect the binding of DNA polymerase to the ssDNA
    The size of the bubble affect the binding
    k_binding comes from this paper

"""


from GIP_rates import replicationRate
from Bio.Seq import Seq
import rxns_CME


def addRepInit(sim,sim_properties):

    n_rxns_init = 0; subsystem = 'initiation'
    # Callled every minute
    
    cell_volume = sim_properties['volume_L'][-1]

    NA = 6.022e23 # Avogadro's

    k_high = 7800*1000/NA/cell_volume # 0.386 molecule^-1 sec^-1
    k_low = 35*1000/NA/cell_volume  # 0.0017 molecule^-1 sec^-1

    # k_on = 100*1000/NA/cell_volume # 0.005 molecule^-1 sec^-1 Number of P_0001/dnaA is above 100
    # k_off = 0.55 #sec^-1

    k_on = 140*1000/NA/cell_volume # 0.005 molecule^-1 sec^-1 Number of P_0001/dnaA is above 100
    k_off = 0.42 #sec^-1

    helicase_removal_rate = 600 #s^-1  SUPER FAST REMOVAL
    
    # k_binding = 2.5e5/NA/cell_volume
    k_replisome_binding = 1E6/NA/cell_volume

    ini_list = []
    ini_counts = []

    # replisomes = ['replisome1', 'replisome2']
    # iniNums_replisomes = [6,12]

    # ini_list.extend(replisomes)
    # ini_counts.extend(iniNums_replisomes)

    # sim.defineSpecies(replisomes)

###########################################################################################################################3
    
    nonOricSpecs = ['High_Affinity_Site','High_Affinity_Bound','High_Affinity_Site_oriC',
                   'High_Affinity_Bound_oriC','Low_Affinity_Site_1','Low_Affinity_Site_2',
                   'Low_Affinity_Bound_1','Low_Affinity_Bound_2','chromosome_C','chromosome_CC']
    
    # High_Affinity_Site are other sites nonrelevant to DNA replication and are competitive to the binding of Site_oriC
    iniNums_nonOricSpec = [16,0,1,0,0,0,0,0,1,1]

    ini_list.extend(nonOricSpecs)
    ini_counts.extend(iniNums_nonOricSpec)

    sim.defineSpecies(nonOricSpecs)

    sim.addReaction(('High_Affinity_Site', 'P_0001'),'High_Affinity_Bound',k_high); n_rxns_init +=1
    sim.addReaction(('High_Affinity_Site_oriC', 'P_0001'),('High_Affinity_Bound_oriC','Low_Affinity_Site_1'),k_high); n_rxns_init +=1
    sim.addReaction(('Low_Affinity_Site_1', 'P_0001'),('Low_Affinity_Bound_1','Low_Affinity_Site_2'),k_low); n_rxns_init += 1
    sim.addReaction(('Low_Affinity_Site_2', 'P_0001'),('Low_Affinity_Bound_2','ssdnaAFila_0'),k_low); n_rxns_init +=1


    unbound_unbind = []
    unbound_unbind_counts = []
    for i in range(31):         #loop adds 30 terms for unbound sites
        unbound = 'ssdnaAFila_' + str(i)
        unbound_unbind.append(unbound)
        unbound_unbind_counts.append(0)
    for k in range(31):
        unbind = 'ssDNA_Unbinding_' + str(k)
        unbound_unbind.append(unbind)
        unbound_unbind_counts.append(0)
    for i in range(20,31):
        unbound_unbind.append('ssdnaAFila_' + str(i)+'_replisome')
        unbound_unbind_counts.append(0)

    unbound_unbind.append('Initiator_C')
    unbound_unbind_counts.append(0)
    unbound_unbind.append('Initiator_CC')
    unbound_unbind_counts.append(0)
    # Marking species
    unbound_unbind.append('RepInitCheck')
    unbound_unbind_counts.append(0)

    ini_list.extend(unbound_unbind)
    ini_counts.extend(unbound_unbind_counts)


    sim.defineSpecies(unbound_unbind)
    

    for i in range (0,30):
        sim.addReaction(('ssdnaAFila_' + str(i),'P_0001'), 'ssdnaAFila_' + str(i+1),k_on); n_rxns_init +=1
        sim.addReaction('ssdnaAFila_' + str(i+1),('ssdnaAFila_' + str(i),'P_0001'),k_off); n_rxns_init +=1

    for i in range(20,31):
        sim.addReaction(('ssdnaAFila_' + str(i), 'P_0044'), ('ssdnaAFila_' + str(i)+'_replisome'), k_replisome_binding); n_rxns_init +=1
        sim.addReaction(('ssdnaAFila_' + str(i)+'_replisome', 'P_0044'), ('ssDNA_Unbinding_' + str(i),'Initiator_C','Initiator_CC','RepInitCheck'), k_replisome_binding); n_rxns_init +=1
        
        # sim.addReaction(('ssdnaAFila_' + str(i), 'replisome1'),('ssDNA_Unbinding_' + str(i),'Initiator_C','Initiator_CC','RepInitCheck'),k_binding); n_rxns_init +=1

    # Add the unbinding reactions for each of the 30 possible unbinding events in filament formation.
    for i in range (1,31):
        sim.addReaction('ssDNA_Unbinding_' + str(i),('ssDNA_Unbinding_' + str(i-1),'P_0001'),helicase_removal_rate); n_rxns_init +=1

    


    
# #################################################################################################################################
        
#     nonOricSpecs2 = ['High_Affinity_Site2','High_Affinity_Bound2','High_Affinity_Site_oriC2',
#                    'High_Affinity_Bound_oriC2','Low_Affinity_Site2_1','Low_Affinity_Site2_2',
#                    'Low_Affinity_Bound2_1','Low_Affinity_Bound2_2']
    
#     sim.defineSpecies(nonOricSpecs2)
#     iniNums_nonOricSpecs2 = [0,0,0,0,0,0,0,0]
    
#     ini_list.extend(nonOricSpecs2)
#     ini_counts.extend(iniNums_nonOricSpecs2)

#     sim.addReaction(('High_Affinity_Site2', 'P_0001'),'High_Affinity_Bound2',k_high); n_rxns_init +=1
#     sim.addReaction(('High_Affinity_Site_oriC2', 'P_0001'),('High_Affinity_Bound_oriC2','Low_Affinity_Site2_1'),k_high); n_rxns_init +=1
#     sim.addReaction(('Low_Affinity_Site2_1', 'P_0001'),('Low_Affinity_Bound2_1','Low_Affinity_Site2_2'),k_low); n_rxns_init +=1
#     sim.addReaction(('Low_Affinity_Site2_2', 'P_0001'),('Low_Affinity_Bound2_2','ssdnaAFila2_0'),k_low); n_rxns_init +=1

#     unbound_unbind2 = []
#     unbound_unbind2_counts = []
#     for i in range(31):         #loop adds 30 terms for different filament length
#         unbound2 = 'ssdnaAFila2_' + str(i)
#         unbound_unbind2.append(unbound2)
#         unbound_unbind2_counts.append(0)
#     for k in range(31):
#         unbind2 = 'ssDNA_Unbinding2_' + str(k)
#         unbound_unbind2.append(unbind2)
#         unbound_unbind2_counts.append(0)

#     unbound_unbind2.append('Initiator_C2')
#     unbound_unbind2_counts.append(0)
#     unbound_unbind2.append('Initiator_CC2')
#     unbound_unbind2_counts.append(0)
#     # Marking species
#     # unbound_unbind.append('RepInitCheck')

#     sim.defineSpecies(unbound_unbind2)
#     ini_list.extend(unbound_unbind2)
#     ini_counts.extend(unbound_unbind2_counts)

#     for i in range (0,30):
#         sim.addReaction(('ssdnaAFila2_' + str(i),'P_0001'), 'ssdnaAFila2_' + str(i+1),k_on); n_rxns_init +=1
#         sim.addReaction('ssdnaAFila2_' + str(i+1),('ssdnaAFila2_' + str(i),'P_0001'),k_off); n_rxns_init +=1
    
#     for i in range(20,31):
#         sim.addReaction(('ssdnaAFila2_' + str(i), 'replisome2'),('ssDNA_Unbinding2_' + str(i),'Initiator_C','Initiator_CC','RepInitCheck'),k_binding); n_rxns_init +=1

#     # Add the unbinding reactions for each of the 30 possible unbinding events in filament formation.
#     for i in range (1,31):
#         sim.addReaction('ssDNA_Unbinding2_' + str(i),('ssDNA_Unbinding2_' + str(i-1),'P_0001'),helicase_removal_rate); n_rxns_init +=1




# #################################################################################################################################
#     nonOricSpecs3 = ['High_Affinity_Site3','High_Affinity_Bound3','High_Affinity_Site_oriC3',
#                    'High_Affinity_Bound_oriC3','Low_Affinity_Site3_1','Low_Affinity_Site3_2',
#                    'Low_Affinity_Bound3_1','Low_Affinity_Bound3_2']
#     sim.defineSpecies(nonOricSpecs3)
#     iniNums_nonOricSpecs3 = [0,0,0,0,0,0,0,0]

#     ini_list.extend(nonOricSpecs3)
#     ini_counts.extend(iniNums_nonOricSpecs3)

#     sim.addReaction(('High_Affinity_Site3', 'P_0001'),'High_Affinity_Bound3',k_high); n_rxns_init +=1
#     sim.addReaction(('High_Affinity_Site_oriC3', 'P_0001'),('High_Affinity_Bound_oriC3','Low_Affinity_Site3_1'),k_high); n_rxns_init +=1
#     sim.addReaction(('Low_Affinity_Site3_1', 'P_0001'),('Low_Affinity_Bound3_1','Low_Affinity_Site3_2'),k_low); n_rxns_init +=1
#     sim.addReaction(('Low_Affinity_Site3_2', 'P_0001'),('Low_Affinity_Bound3_2','ssdnaAFila3_0'),k_low); n_rxns_init +=1

#     unbound_unbind3 = []
#     unbound_unbind3_counts = []

#     for i in range(31):         #loop adds 30 terms for unbound sites
#         unbound3 = 'ssdnaAFila3_' + str(i)
#         unbound_unbind3.append(unbound3)
#         unbound_unbind3_counts.append(0)

#     for k in range(31):
#         unbind3 = 'ssDNA_Unbinding3_' + str(k)
#         unbound_unbind3.append(unbind3)
#         unbound_unbind3_counts.append(0)

#     unbound_unbind3.append('Initiator_C3')
#     unbound_unbind3_counts.append(0)

#     unbound_unbind3.append('Initiator_CC3')
#     unbound_unbind3_counts.append(0)
#     # Marking species
#     # unbound_unbind.append('RepInitCheck')

#     sim.defineSpecies(unbound_unbind3)
#     ini_list.extend(unbound_unbind3)
#     ini_counts.extend(unbound_unbind3_counts)

#     for i in range (0,30):
#         sim.addReaction(('ssdnaAFila3_' + str(i),'P_0001'), 'ssdnaAFila3_' + str(i+1),k_on); n_rxns_init +=1
#         sim.addReaction('ssdnaAFila3_' + str(i+1),('ssdnaAFila3_' + str(i),'P_0001'),k_off); n_rxns_init +=1
    
#     for i in range(20,31):
#         sim.addReaction(('ssdnaAFila3_' + str(i), 'replisome2'),('ssDNA_Unbinding3_'+str(i),'Initiator_C','Initiator_CC','RepInitCheck'),k_binding); n_rxns_init +=1

#     # Add the unbinding reactions for each of the 30 possible unbinding events in filament formation.
#     for i in range (1,31):
#         sim.addReaction('ssDNA_Unbinding3_' + str(i),('ssDNA_Unbinding3_' + str(i-1),'P_0001'),helicase_removal_rate); n_rxns_init +=1

#     print('Total {0} reaction added into replication initiation'.format(n_rxns_init))

    for i in range(n_rxns_init):
        rxns_CME.addReactionToMap(sim_properties, subsystem)

    return ini_list, ini_counts



def addReplication(sim,sim_properties,gbfile):   
    """
    

    Description: 3 pseudo genes are not included in replication
    """ 
    subsystem = 'replication'

    dna = gbfile

    # sim.defineSpecies(['G_0051', 'G_0602', 'G_0546'])

    gene_list = []
    # 493 elements in gene_list

    for i in range(len(dna.features)):
        if ('product' in dna.features[i].qualifiers.keys()):
            #print(i) # This first statement works
            #print(dna.features[i].qualifiers['product'])
            if dna.features[i].qualifiers['product'][0]:# Figure out how to sort out for ribosomal operons?
                #print(dna.features[i].qualifiers['product'])
                gene_list.append(i)

    # print(gene_list)
    # [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, ...]
                
    # rep_list to store the names of all intergenic species and 
    # will be called at the end of this function to give the numbers of each intergenic species
    rep_list = []
    rep_counts = []

    C_bp = 0
    CC_bp = 0

    position = 0
    
    CC_genes = []

    for gene in gene_list:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real

        if start < len(dna.seq)/2:
            geneSeq = Seq(str(dna.seq[position:end]))#, generic_dna)
    #         print(locusTag)

            if start == 0:
                #print(locusTag)

                n_tot, k_gene_rep = replicationRate(sim_properties, locusNum)
                intergenic_region = locusTag+'_inter'

                geneID = 'G_' + locusNum
                GeneProduced = 'Produced_G_' +locusNum

                RepProd = [geneID,GeneProduced,intergenic_region]

                sim.defineSpecies(RepProd)

                sim.addReaction(('chromosome_C','Initiator_C'), tuple(RepProd), k_gene_rep)
                rxns_CME.addReactionToMap(sim_properties, subsystem, locusNum=locusNum )

                rep_list.extend(RepProd)
                rep_counts.extend([1,0,0])
                
                position  = end
                C_bp += n_tot

            else:

                #print(locusTag)
                #print(position,end)
                n_tot, k_gene_rep = replicationRate(sim_properties, locusNum)
            
                intergenic_region = locusTag+'_inter' 
                geneID = 'G_' + locusNum
                GeneProduced = 'Produced_G_' +locusNum
                RepProd = [geneID,GeneProduced,intergenic_region]

                sim.defineSpecies(RepProd)


                sim.addReaction(rep_list[-1], tuple(RepProd), k_gene_rep)
                rxns_CME.addReactionToMap(sim_properties, subsystem, locusNum=locusNum )
                
                rep_list.extend(RepProd)
                rep_counts.extend([1,0,0])

                position = end
                C_bp += n_tot


    for gene in gene_list:
        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real

        if start > len(dna.seq)/2:
    
            CC_genes.append(gene)

#     print(CC_genes)
    CC_genes.reverse()
#     print(CC_genes)
    # The very end of the chromosome

    position = 543086

    for gene in CC_genes:

        locusTag = dna.features[gene].qualifiers['locus_tag'][0]
        locusNum = locusTag.split('_')[1] 
        start =  dna.features[gene].location.start.real
        end  = dna.features[gene].location.end.real

        geneSeq = Seq(str(dna.seq[start:position]))#, generic_dna)
    #         print(locusTag)

        if end == 543086:
            #print(locusTag)
            n_tot, k_gene_rep = replicationRate(sim_properties, locusNum)

            intergenic_region = locusTag+'_inter'

            geneID = 'G_' + locusNum

            GeneProduced = 'Produced_G_' +locusNum

            RepProd = [geneID,GeneProduced,intergenic_region]

            sim.defineSpecies(RepProd)

            sim.addReaction(('chromosome_CC','Initiator_CC'), tuple(RepProd), k_gene_rep)
            rxns_CME.addReactionToMap(sim_properties, subsystem, locusNum=locusNum )

            rep_list.extend(RepProd)
            rep_counts.extend([1,0,0])   

            position  = start
            CC_bp += n_tot

        # This if elif else conditionals are tricky because 
        elif dna.features[gene].qualifiers['locus_tag'][0] == 'JCVISYN3A_0421':

            # print('End of Replication')

            n_tot, k_gene_rep = replicationRate(sim_properties, locusNum)

            intergenic_region = locusTag+'_inter' 

            geneID = 'G_' + locusNum

            GeneProduced = 'Produced_G_' +locusNum
            RepProd = [geneID,GeneProduced,intergenic_region]

            sim.defineSpecies(RepProd)
            
            gene_rep_end_products = ['chromosome_C','chromosome_CC','chromosome_C','chromosome_CC', 'P_0044', 'P_0044'] # Return P_0044 after replication

            gene_rep_end_products.extend(RepProd)

            sim.addReaction(rep_list[-1], tuple(gene_rep_end_products), k_gene_rep)
            rxns_CME.addReactionToMap(sim_properties, subsystem, locusNum=locusNum )

            position = start
            CC_bp += n_tot
            rep_list.extend([geneID, GeneProduced, intergenic_region])
            rep_counts.extend([1,0,0])

        else:

            #print(locusTag)
            #print(position,end)
            n_tot, k_gene_rep = replicationRate(sim_properties, locusNum)

            intergenic_region = locusTag+'_inter' 

            geneID = 'G_' + locusNum
            GeneProduced = 'Produced_G_' +locusNum
            RepProd = [geneID,GeneProduced,intergenic_region]

            sim.defineSpecies(RepProd)

            sim.addReaction(rep_list[-1], tuple(RepProd), k_gene_rep)
            rxns_CME.addReactionToMap(sim_properties, subsystem, locusNum=locusNum)

            position = start
            CC_bp += n_tot
            rep_list.extend(RepProd)
            rep_counts.extend([1,0,0])




    return rep_list, rep_counts