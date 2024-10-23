import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, PeriodicTable, Atom, rdMolDescriptors

from mordred import Calculator, RingCount, BondCount, AdjacencyMatrix
from mordred._atomic_property import *

from descriptors.utils.chem import *

from math import sqrt, pi, log

from descriptors.utils.chem import GetPrincipleQuantumNumber


def atom_properties(mol, **kwargs): # revised in v3.0.
    # Ensure the property is implemented
    try:
        mol = Chem.AddHs(mol) # H added to fix error in calculation. Version v3.0 wotan
        property_fn = ATOM_PROPERTIES[kwargs['prop']]
        values = [property_fn(atom) for atom in mol.GetAtoms()]


        # Different behaviour depending if they are atom or molecule properties
        # atom_property = kwargs['prop'] not in ['bo', 'kh']
        # if atom_property:
        #     values = [property_fn(atom) for atom in mol.GetAtoms()]
        #
        # else:
        #     values = property_fn(mol)

    except KeyError:
        raise NotImplementedError


    # Scale with respect to cabron atom if required
    if kwargs['scale']:
        values = [val / property_fn(Chem.Atom(6)) for val in values]



    # Apply corresponding transformations
    if kwargs['type'] == 'sum':
        return sum(values)

    elif kwargs['type'] == 'mean':
        return sum(values) / len(values)

    else:
        raise NotImplementedError


def mol_properties(mol, **kwargs):
    # Ensure the property is implemented
    try:
        return MOL_PROPERTIES[kwargs['prop']](mol, **kwargs)

    except KeyError:
        raise NotImplementedError



##############################
###  Molecule descriptors  ###
##############################

def MolecularWeight(mol, **kwargs): # revised in v3.0.
    '''Molecular weight of the molecule (average or not)
    '''
    if kwargs['type'] == 'mean':

        mol = Chem.AddHs(mol)
        # return Descriptors.MolWt(mol) / mol.GetNumAtoms() # changed to be comparable to PaDEL.  Version v3.0 wotan
        return Descriptors.ExactMolWt(mol) / mol.GetNumAtoms() #https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html

    elif kwargs['type'] == '':
        mol = Chem.AddHs(mol)
        # return Descriptors.MolWt(mol) # changed to be comparable to PaDEL. Version v3.0 wotan
        return Descriptors.ExactMolWt(mol) #https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html


# def KierHallAtom(atom): # not implemented as individual descriptor. used to KierHallMolecule
#     '''Kier-Hall electrotopological state for a single atom
#     '''
#     atomic_num = atom.GetAtomicNum()
#
#     if atomic_num == 1:
#         return -0.2
#     else:
#         # TODO: Check with Scikit-Chem
#         # NOTES: when revised, in rdkit it is described as:
#         #         atom.GetImplicitValence : Returns the explicit valence of the atom.
#         #         atom.GetExplicitValence : Returns the number of implicit Hs on the atom.
#         #        Changed in order to make the same calculation as in Maestro
#         return float(atom.GetImplicitValence() - atom.GetExplicitValence()) / \
#                (GetPrincipleQuantumNumber(atomic_num) ** 2)




def KierHallAtom(atom): # revised in v3.0
    '''Kier-Hall electrotopological state for a single atom
    '''
    atomic_num = atom.GetAtomicNum()

    delta_vi = ValenceVertexDegree2(atom)
    delta_i = atom.GetDegree() # changed accordingly to Maestro Schödinger
    L_i = GetPrincipleQuantumNumber(atomic_num)

    if atomic_num == 1:
        return -0.2
    else:
        # TODO: Check with Scikit-Chem

        return float(delta_vi - delta_i)/(L_i ** 2)


def KierHallMolecule(mol, **kwargs): # disapeared after Dragon 5.5
    '''Kier-Hall electrotopological states for a molecule

    NOTE: Returns Sum similar to Maestro-Schödinger but very
        dissimilar to Dragon
    '''
    mol = Chem.RemoveHs(mol) # added in wotan v3.0

    if kwargs['type'] == 'sum':
        return sum([KierHallAtom(atom) for atom in mol.GetAtoms()])

    elif kwargs['type'] == 'mean':
        return sum([KierHallAtom(atom) for atom in mol.GetAtoms()]) \
               / mol.GetNumAtoms()


def ConventionalBondOrders(mol, **kwargs):
    '''The sum of conventional bond orders (SCBO) is obtained by adding up the
    conventional bond order of all the edges in the H-depleted molecular graph
    representing a molecule. Conventional bond orders are: 1 for simple bond,
    2 for double bond, 3 for triple bond and 1.5 for aromatic bond.
    '''
    '''
        NOTES to v3.0:
        exist a little discrepancy with dragon, dragon7 and Vega due to its different
        criteria on single bonds (see Bonds function).
    '''
    mol = Chem.RemoveHs(mol)

    if kwargs['type'] == 'sum':
        return sum([bond.GetBondTypeAsDouble() for bond in mol.GetBonds()])

    elif kwargs['type'] == 'mean':
        return sum([bond.GetBondTypeAsDouble() for bond in mol.GetBonds()]) \
               / mol.GetNumAtoms()


def TotalAtomNumber(mol, **kwargs):
    '''Compute the number of atoms including (kwargs['hydrogens'] = true)
       or not including hydrogens (kwargs['hydrogens'] = false)
    '''

    if kwargs['hydrogens']:
        mol = Chem.AddHs(mol)

    return mol.GetNumAtoms()

def AtomNumber(mol, **kwargs):
    '''Compute the number of atoms with kwargs[atom]
    '''

    mol = Chem.AddHs(mol)

    smart = Chem.MolFromSmarts(ATOMS[kwargs['atom']])
    matches = len(mol.GetSubstructMatches(smart))

    return matches


def Bonds(mol, **kwargs):
    '''Calculates the number of bonds (non H) and the multiple bond Number

        kwargs[type] = 'any' - all bonds
        kwargs[type] = 'heavy' - all non-H bonds
        kwargs[type] = 'multiple' - all multiple bonds
        kwargs[type] = 'single'

        kwargs[type] = 'double'

        kwargs[type] = 'triple'
        kwargs[type] = 'aromatic'

        all them are calculated in non-ketulized structure.
        Are different from those calculated in dragon, padel
    '''

    calc = Calculator(
        BondCount.BondCount(kwargs['type'], False)
    )

    desc = calc.pandas([mol], nproc=1, quiet=True)

    return desc.iloc[0, 0]


def ARR(mol, **kwargs):
    '''Calculates the ratio of aromatic bonds over total number of non-H (heavy) GetNumBonds
    '''
    ''' v3.0 has changed to avoid smarts substructure search
    '''
    calc_arom_bond = Calculator(
        BondCount.BondCount('aromatic', False)
    )

    desc_arom_bond = calc_arom_bond.pandas([mol], nproc=1, quiet=True)

    arom_bonds = desc_arom_bond.iloc[0, 0]

    calc_heavy = Calculator(
        BondCount.BondCount('heavy', False)
    )

    desc_heavy = calc_heavy.pandas([mol], nproc=1, quiet=True)

    heavy_bonds = desc_heavy.iloc[0, 0]

    # bond_number = mol.GetNumBonds()
    #
    # arom_bond_smart = Chem.MolFromSmarts('*:*')
    # arom_bonds = len(mol.GetSubstructMatches(arom_bond_smart))

    if heavy_bonds != 0:
        return arom_bonds / heavy_bonds
    else:
        return 0


def Rings(mol, **kwargs):
    '''Calculates the number of rings in the molecule
    '''
    return rdMolDescriptors.CalcNumRings(mol)


def RBN(mol, **kwargs):
    '''Calculates the number of rotable bonds (and fraction) in the molecule
    '''
    '''
    is equivalent to nRotB and nRtotBFrac in PaDEL(that excludes terminal
    bonds), and same names in Dragon,  but gives different results
    '''

    # added parameter strict = False in v3.0 makes results equal to cdk and Padel.

    rotable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol, strict = False)

    mol = Chem.RemoveHs(mol)
    bond_number = mol.GetNumBonds()

    if kwargs['type'] == 'fraction':
        if rotable_bonds == 0:
            return 0
        else:
            return rotable_bonds/bond_number

    else:
        return rotable_bonds


def nR(mol, **kwargs):
    '''Number of [length]-membered rings

    Keyword arguments:
        length -- length of the cycle
        greater -- count length or greater rings
        fused -- count fused rings
        aromatic -- count aromatic (true), non-aromatic (false) or both (none)
        hetero -- count hetero (true), carbon (false) or both rings (None)
    '''

    for argument,value in kwargs.items():
        if value == 'None':
            kwargs[argument] = None

    calc = Calculator(
        RingCount.RingCount(kwargs['length'], kwargs['greater'],
                            kwargs['fused'], kwargs['aromatic'],
                            kwargs['hetero'],)
    )

    desc = calc.pandas([mol], nproc=1, quiet=True)

    return desc.iloc[0, 0]


# TODO: Add function to compute all the descriptors, dictionaries with
# descriptor name and kwargs maybe? MAAAAAAAAAAAYBE


ATOM_PROPERTIES = {
    'vdw': get_vdw_volume,
    'eneg': get_sanderson_en,
    'polar': get_polarizability,
}

MOL_PROPERTIES = {
    'mw': MolecularWeight,
    'kh': KierHallMolecule,
    'bo': ConventionalBondOrders,
    'tan': TotalAtomNumber,
    'an': AtomNumber,
    'nb': Bonds,
    'ar': ARR,
    'rg': Rings,
    'rb': RBN,
    'nr': nR,
}

ATOMS = {
    'H': '[#1]',
    'C': '[#6]',
    'N': '[#7]',
    'O': '[#8]',
    'P': '[#15]',
    'S': '[#16]',
    'F': '[#9]',
    'CL': '[#17]',
    'BR': '[#35]',
    'I': '[#53]',
    'B': '[#5]',
}


def Constitutional_all(mol):

    desc_df = pd.DataFrame()

    types = ['sum', 'mean']
    for key in ATOM_PROPERTIES.keys():
        for type in types:
            kwargs = {'prop': key, 'scale': True, 'type': type}
            desc_name = '{}{}'.format(type[0].upper(), key[0])
            desc_value = atom_properties(mol, **kwargs)
            desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

            desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)


#    MolecularWeight
    for type in ['', 'mean']:
        kwargs = {'prop': 'mw', 'type': type}
        desc_name = 'AMW' if type else 'MW'
        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
#
#
#    # KierHallMolecule
    for type in ['sum', 'mean']:
        kwargs = {'prop': 'kh', 'type': type}
        desc_name = '{}s'.format(type[0].upper())
        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
#
#     #
#     # # NOTE: this one is yielding empty dataframes
    # # TotalAtomNumber (With and w/o hydrogens)
    for type in ['', 'hydrogens']:
        kwargs = {'prop': 'tan', 'hydrogens': type}
        desc_name = 'nAT' if type else 'nSK'
        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
#
    # Atomnumber
    for atom in ATOMS.keys():
        kwargs = {'prop': 'an', 'atom': atom}
        desc_name = 'n{}'.format(atom)
        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
#
#
#     # Bonds
    bond_types = \
        ['any', 'heavy', 'multiple', 'single', 'double', 'triple', 'aromatic']
    for type in bond_types:
        kwargs = {'prop': 'nb', 'type': type}

        if   type == 'any':     desc_name = 'nBT'
        elif type == 'heavy':   desc_name = 'nBO'
        else:                   desc_name = 'n' + type.capitalize()

        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
#
#
    # # RBN
    for type in ['', 'fraction']:
        kwargs = {'prop': 'rb', 'type': type}
        desc_name = 'RBF' if type else 'RBN'
        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
#
#     # Check nR names to adapt (some of them) to Dragon
#     # nR
    for length in range(3, 13):
        for greater in [True, False]:
            for fused in [True, False]:
                for aromatic in [True, False, None]:
                    for hetero in [True, False, None]:
                        kwargs = {
                            'prop': 'nr',
                            'length': length,
                            'greater': greater,
                            'fused': fused,
                            'aromatic': aromatic,
                            'hetero': hetero,
                        }
                        desc_name = 'nR_{}_{}_{}_{}_{}'.format(length, greater, fused, aromatic, hetero)
                        desc_value = mol_properties(mol, **kwargs)
                        desc_value = \
                            pd.DataFrame.from_dict({desc_name: [desc_value]})

                        desc_df = \
                            pd.concat([desc_df, desc_value], axis=1, sort=False)
#
#     #
#     # ConventionalBondOrders
    for type in ['sum', 'mean']:
        kwargs = {'prop': 'bo', 'type': type}
        desc_name = '{}CBO'.format(type[0].upper())
        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
#     #
#
    # Properties without keyword arguments
    names_dict = {'ar': 'ARR', 'rg': 'nCIC'}
    for prop in ['ar', 'rg']:
        kwargs = {'prop': prop}
        desc_name = names_dict[prop]
        desc_value = mol_properties(mol, **kwargs)
        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

    return desc_df
