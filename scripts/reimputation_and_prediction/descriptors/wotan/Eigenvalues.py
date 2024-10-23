

from mordred import Calculator
from mordred import BaryszMatrix

from rdkit import Chem
from rdkit.Chem import rdmolops

import numpy as np
import pandas as pd


def LP1(mol, **kwargs):
    '''
    leading eigenvalue from adjacency matrix (Lovasz-Pelikan index)
    called LP1 in dragon <5.5 and equivalent to SpMAX_A in dragon7'''
    adj_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol)
    eigenvals = np.linalg.eigvalsh(adj_matrix)

    return eigenvals[-1]


def Barysz(mol, **kwargs):
    '''
    https://mordred-descriptor.github.io/documentation/master/matrix.html#matrix-aggregating-methods
    SpAbs : graph energy from Barysz matrix
    SpMax : normalized leading eigenvalue from Barysz matrix
    SpDiam : spectral diameter from Barysz matrix
    SpAD : spectral absolute deviation from Barysz matrix
    SpMAD : spectral mean absolute deviation from Barysz matrix
    LogEE : Estrada-like index (log function) from Barysz matrix
    SM1 - SM3 : spectral moment of order 1-3
    VE1 : coefficient sum of the last eigenvector from Barysz matrix
    VE2 : average coefficient of the last eigenvector from Barysz matrix
    VE3 : logarithmic coefficient sum of the last eigenvector from Barysz matrix
    VR1 : Randic-like eigenvector-based index from Barysz matrix
    VR2 : normalized Randic-like eigenvector-based index from Barysz matrix
    VR3 : logarithmic Randic-like eigenvector-based index from Barysz matrix

    all can be weigthened by:
    https://mordred-descriptor.github.io/documentation/master/atomic_prop.html#atomic-properties
        Z : atomic number
        m : mass
        v : van der Waals volume
        se : sanderson EN
        pe : pauling EN
        are : allred-rocow EN
        p : polarizability
        i : ionization potential
    '''
    calc = Calculator(
    BaryszMatrix.BaryszMatrix(kwargs['prop'], kwargs['type'])
    )

    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        return(desc.iloc[0, 0])

    else: return np.nan


def Eigenvalues_all(mol): # all checked in v3.0. concordant with mordred online


    properties = ['Z', 'm', 'v', 'se', 'pe', 'are', 'p', 'i']
    types = ['SpAbs', 'SpMax', 'SpDiam', 'SpAD', 'SpMAD', 'LogEE', 'SM1', 'VE1', 'VE2', 'VE3', 'VR1', 'VR2', 'VR3']

    desc_df = pd.DataFrame()

    desc_name = 'LP1'
    desc_value = LP1(mol)
    desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

    desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

    for type in types:
        for prop in properties:
            kwargs = {'prop': prop, 'type': type}

            desc_name = '{}_{}'.format(type, prop) # in mordred: '{}_Dz{}'
            desc_value = Barysz(mol, **kwargs)
            desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

            desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

    return desc_df
