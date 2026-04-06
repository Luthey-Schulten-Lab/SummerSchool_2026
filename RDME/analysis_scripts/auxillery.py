from math import floor
from math import log10
import numpy as np
def round_sig(x, sig=2):
    negative = False
    if x < 0:
        negative = True
    x = abs(x)
    if negative:
        return -1*round(x, sig-int(floor(log10(abs(x))))-1)
    elif x==0.0:
        return 0.0
    else:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    
    
def getAvgRxnTraceCounts(rxnID, w):
    allFluxes = w.get_rxn_trace(rxnID)
    allFluxes = allFluxes.T

    vol = ['Volume']
    
    vols = w.get_species_traces(vol)
    vols = vols[0].T/1000000*1e-15
#     print(vols.shape)
    
    allCountFluxes = np.multiply(allFluxes,vols)*6.02e23/1000

    
    return allCountFluxes