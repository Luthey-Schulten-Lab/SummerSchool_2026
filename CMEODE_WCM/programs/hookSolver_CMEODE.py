"""
Author: Enguang Fu

Date: March 2024

user defined solver where hookSimulation performs every hookInterval
"""

import communicate
import hook_CMEODE
import GIP_rates

# import pyLM.units to make lm.GillespieDSolver readabe by python
from pyLM.units import *

class MyOwnSolver(lm.GillespieDSolver): 
    """
    Description:
    lm.GillespieDSovler is the parent class or base class and MyOwnSolver is the subclass, thus MyOwnsolver inheritate the getSpeciesCountView() method
    
    """
    def initializeSolver(self, sim, CME_count_array, sim_properties,genome3A):
        """
        Input:
        sim, an instance of Class CME.CMESilumation and its attributes contains the simulation information
        
        """
        self.sim = sim
        self.CME_count_array = CME_count_array
        self.sim_properties = sim_properties
        self.genome3A = genome3A

    def update_countarray(self):
        self.CME_count_array.getSpeciesCountView(self) # Here the self points to MyOwnSolver, a subclass of lm.GillespieDSolver

    def update_rateconstants(self):
        """
        Description: Update the rate constants of replication, trascription and translation reactions
        """
        rxns_map = self.sim_properties['rxns_map']
        
        for rxns_id, i_rxns in rxns_map.items():
            i_rxns = i_rxns # In C++, reaction starts from 1

            # if i_rxns == 0:
            #      print(f"Reaction {i_rxns} {self.getReactionRateConstantsView(i_rxns)}")

            # print out 
            if rxns_id == 'rep_0001':
                    print('Gene 0001 Replication rate {0} before change'.format(self.getReactionRateConstantsView(i_rxns)))
            elif rxns_id == 'trsc_poly_0001':
                    print('mRNA 0001 Transcription rate {0} before change'.format(self.getReactionRateConstantsView(i_rxns)))
            elif rxns_id == 'tran_poly_0001':
                    print('protein 0001 Translation rate {0} before change'.format(self.getReactionRateConstantsView(i_rxns)))
            else:
                 None

            rxns_prefixes = ['rep', 'trsc_poly', 'tran_poly', 'deg_depoly']

            # update CME rates per second after ODE
            if rxns_id.startswith('rep'):
                rate_constant = self.getReactionRateConstantsView(i_rxns) # [rate]
                locusNum = rxns_id.split('_')[-1]
                rate_constant[0] = GIP_rates.replicationRate(self.sim_properties, locusNum)[1]
                if rxns_id == 'rep_0001':
                    print('Gene 0001 Replication rate {0} after change'.format(self.getReactionRateConstantsView(i_rxns)))

            elif rxns_id.startswith('trsc_poly'):
                rate_constant = self.getReactionRateConstantsView(i_rxns) # [rate]
                locusNum = rxns_id.split('_')[-1]
                rate_constant[0] = GIP_rates.TranscriptionRate(self.sim_properties, locusNum)
                if rxns_id == 'trsc_poly_0001':
                    print('mRNA 0001 Transcription rate {0} after change'.format(self.getReactionRateConstantsView(i_rxns)))

            elif rxns_id.startswith('tran_poly'):
                rate_constant = self.getReactionRateConstantsView(i_rxns) # [rate]
                locusNum = rxns_id.split('_')[-1]
                full_aasequence = self.sim_properties['genome']['JCVISYN3A_'+locusNum]['AAsequence']

                rate_constant[0] = GIP_rates.TranslationRate(self.sim_properties, locusNum, full_aasequence)
                if rxns_id == 'tran_poly_0001':
                    print('protein 0001 Translation rate {0} after change'.format(self.getReactionRateConstantsView(i_rxns)))

            # elif rxns_id.startswith('translocation'):
            #     rate_constant = self.getReactionRateConstantsView(i_rxns) # [rate]
            #     locusNum = rxns_id.split('_')[-1]
            #     rate_constant[0] = GIP_rates.TranslocationRate(self.sim_properties, locusNum)
            
            elif rxns_id.startswith('deg_depoly'):
                rate_constant = self.getReactionRateConstantsView(i_rxns) # [rate]
                locusNum = rxns_id.split('_')[-1]
                rate_constant[0] = GIP_rates.mrnaDegradationRate(self.sim_properties, locusNum)

            else:
                None
            

    def hookSimulation(self, time): # hookSimulation method here will override the ones in the C++ code
        
        restartNum = self.sim_properties['restartNum'][-1]
        restartInterval = self.sim_properties['restartInterval']
      
        print('*******************************************************************')
        print(f'HookSimulation is called at {restartNum*restartInterval+time:.3f} second')
        
        self.update_countarray()

        # using self.CMECounts to update CME counts and write them into sim_properties which will be used in ODE
        communicate.updateCMEcountsHook(self.sim, self.CME_count_array, self.sim_properties)
        
        # Communicate Cost; Do ODE; Update SA and volume
        hook_CMEODE.hook_CMEODE(self.sim_properties, self.genome3A)
        
        # update counts of certain ODE species that also in CME to self.CMECounts
        communicate.updateODEtoCME(self.sim, self.CME_count_array, self.sim_properties)

        # update CME rates per communication step after ODE
        self.update_rateconstants()

        return 1
        


        


            