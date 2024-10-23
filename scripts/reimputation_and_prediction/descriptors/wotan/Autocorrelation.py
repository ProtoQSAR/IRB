
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from mordred import Calculator
from mordred import Autocorrelation

from math import log
import numpy as np
import pandas as pd

''''
    Different autocorrelation values weigthed by different properties:
    https://mordred-descriptor.github.io/documentation/master/atomic_prop.html#atomic-properties
            c : charge
            dv : valence electrons
            d : sigma electrons
            s : intrinsic state
            Z : atomic number
            m : mass
            v : van der Waals volume
            se : sanderson EN
            pe : pauling EN
            are : allred-rocow EN
            p : polarizability
            i : ionization potential


'''
def ATS(mol, **kwargs):

    '''Leverage-weighted autocorrelation of lag [order] / weighted by [prop]

    Keyword arguments:
        order -- lag of the autocorrelation
        prop -- property to weight on (usually, mass -- 'm')

    '''

    # The result cannot be checked because dragon has changed the calculation
    # of this descriptor in the new versions --> checkef in v3.0: can be
    # checked in d7 from https://ochem.eu


    calc = Calculator(
        Autocorrelation.ATS(kwargs['order'], kwargs['prop'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        result = desc.iloc[0, 0]

    else:
        result = np.nan

    return result


def AATS(mol, **kwargs):
    '''Leverage-weighted autocorrelation of lag [order] / weighted by [prop]

    Keyword arguments:
        order -- lag of the autocorrelation
        prop -- property to weight on (usually, mass -- 'm')
    '''


    calc = Calculator(
        Autocorrelation.AATS(kwargs['order'], kwargs['prop'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        result = desc.iloc[0, 0]

    else:
        result = np.nan

    return result


def ATSC(mol, **kwargs):
    '''Leverage-weighted autocorrelation of lag [order] / weighted by [prop]

    Keyword arguments:
        order -- lag of the autocorrelation
        prop -- property to weight on (usually, mass -- 'm')
    '''

    calc = Calculator(
        Autocorrelation.ATSC(kwargs['order'], kwargs['prop'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        result = desc.iloc[0, 0]

    else:
        result = np.nan

    return result


def AATSC(mol, **kwargs):
    '''Leverage-weighted autocorrelation of lag [order] / weighted by [prop]

    Keyword arguments:
        order -- lag of the autocorrelation
        prop -- property to weight on (usually, mass -- 'm')
    '''

    calc = Calculator(
        Autocorrelation.AATSC(kwargs['order'], kwargs['prop'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        result = desc.iloc[0, 0]

    else:
        result = np.nan

    return result


def MATS(mol, **kwargs):
    '''Leverage-weighted autocorrelation of lag [order] / weighted by [prop]

    Keyword arguments:
        order -- lag of the autocorrelation
        prop -- property to weight on (usually, mass -- 'm')
    '''

    calc = Calculator(
        Autocorrelation.MATS(kwargs['order'], kwargs['prop'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        result = desc.iloc[0, 0]

    else:
        result = np.nan

    return result


def GATS(mol, **kwargs):
    '''Leverage-weighted autocorrelation of lag [order] / weighted by [prop]

    Keyword arguments:
        order -- lag of the autocorrelation
        prop -- property to weight on (usually, mass -- 'm')
    '''

    calc = Calculator(
        Autocorrelation.GATS(kwargs['order'], kwargs['prop'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        result = desc.iloc[0, 0]

    else:
        result = np.nan

    return result


def Autocorrelation_all(mol):
    '''
    Some descriptors does not exist really in mordred, (e.g. MATS0c), but it will
    return constant values (np.values) and will never be selected
    '''
    autocorr_descriptors = {
        'ATS': ATS,
        'AATS': AATS,
        'ATSC': ATSC,
        'AATSC': AATSC,
        'MATS': MATS,
        'GATS': GATS,
    }

    desc_df = pd.DataFrame()

    for desc in autocorr_descriptors.keys():
        if desc == 'ATS' or desc == 'AATS':
            properties = ['dv', 'd', 's', 'Z', 'm', 'v', 'se', 'pe', 'are', 'p', 'i']
        else:
            properties = ['c', 'dv', 'd', 's', 'Z', 'm', 'v', 'se', 'pe', 'are', 'p', 'i']
        if desc == 'MATS' or desc == 'GATS':
            order_values = [i for i in range(1, 9)]
        else:
            order_values = [i for i in range(0, 9)]
        for order in order_values:
            for prop in properties:
                kwargs = {'order': order, 'prop': prop}
                try:
                    desc_name = '{}{}{}'.format(desc, order, prop)
                    desc_value = autocorr_descriptors[desc](mol, **kwargs)
                    desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

                    desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

                except: continue

    return desc_df
