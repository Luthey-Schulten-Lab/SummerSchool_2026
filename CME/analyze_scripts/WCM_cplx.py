"""
Functions used in complex formation analysis
"""

import pandas as pd
import numpy as np


def get_cplx_dict(cplx_path):
    """
    Return a list of complexes involved in the simulation
    """
    
    cplx_df = pd.read_csv(cplx_path, comment='#')

    cplx_dict = {}

    for index, row in cplx_df.iterrows():
        name = row['Name']
        cplx_dict[name] = {}
        locusNums = row['Genes'].split(';')
        cplx_dict[name]['locusNums'] = locusNums
        cplx_dict[name]['init_count'] = row['Count']


    return cplx_dict

