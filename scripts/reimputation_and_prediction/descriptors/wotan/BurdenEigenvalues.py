
from mordred import Calculator, BCUT

import numpy as np
import pandas as pd

def BELpk(mol, **kwargs):
    '''eigenvalue n. k of Burden matrix weighted by property
                    (c, dv d, s, Z, m, v se, pe, are, p i, )
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

## Mordred direct timplementation

    calc = Calculator(
        BCUT.BCUT(kwargs['prop'], kwargs['order'])
    )

    desc = calc.pandas([mol], nproc=1, quiet=True)

    if type(desc.iloc[0, 0]) is np.float64:
        return(desc.iloc[0, 0].round(2))

    else: return np.nan

def BurdenEigenvalues_all(mol):

    order_values = [0, -1]
    properties = ['c', 'dv', 'd', 's', 'Z', 'm', 'v', 'se', 'pe', 'are', 'p', 'i']

    desc_df = pd.DataFrame()

    for prop in properties:
        for order in order_values:
            kwargs = {'prop': prop, 'order': order}

            desc_name = 'BEL{}{}'.format(prop, order)
            desc_value = BELpk(mol, **kwargs)
            desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

            desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

    return desc_df
