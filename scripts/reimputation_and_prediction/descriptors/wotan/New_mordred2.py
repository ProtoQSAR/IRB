"""
Created on Wed Jun 10 12:32:36 2020

@author: Pablo Aparicio
"""

from mordred import Calculator, descriptors
from mordred import MoeType, CPSA, EState, SLogP
from rdkit import Chem
from rdkit.Chem import AllChem, AddHs

import numpy as np
import pandas as pd


######### MoeType descriptors #########

def LabuteASA(mol, **kwargs):


    calc = Calculator(
        MoeType.LabuteASA()
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    return desc.iloc[0,0]


def PEOE_VSA(mol, **kwargs):


    calc = Calculator(
        MoeType.PEOE_VSA(kwargs['k'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    return desc.iloc[0,0]


def SMR_VSA(mol, **kwargs):


    calc = Calculator(
        MoeType.SMR_VSA(kwargs['k'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def SLogP_VSA(mol, **kwargs):


    calc = Calculator(
        MoeType.SlogP_VSA(kwargs['k'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def EState_VSA(mol, **kwargs):


    calc = Calculator(
        MoeType.EState_VSA(kwargs['k'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def VSA_EState(mol, **kwargs):


    calc = Calculator(
        MoeType.VSA_EState(kwargs['k'])
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]



######### CPSA descriptors #########

def PNSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.PNSA(kwargs['version'])
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.PNSA(kwargs['version'])
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]

def PPSA(mol, **kwargs):
    
    # print(Chem.MolToSmiles(mol)) #to debug

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)


    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.PPSA(kwargs['version'])
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.PPSA(kwargs['version'])
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)



    return desc.iloc[0,0]


def DPSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.DPSA(kwargs['version'])
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.DPSA(kwargs['version'])
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)

    return desc.iloc[0,0]


def FNSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.FNSA(kwargs['version'])
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.FNSA(kwargs['version'])
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)



    return desc.iloc[0,0]

def FPSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.FPSA(kwargs['version'])
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.FPSA(kwargs['version'])
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def WNSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.WNSA(kwargs['version'])
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.WNSA(kwargs['version'])
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)



    return desc.iloc[0,0]


def WPSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.WPSA(kwargs['version'])
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.WPSA(kwargs['version'])
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def RNCG(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.RNCG()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.RNCG()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def RPCG(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.RPCG()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.RPCG()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def RNCS(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.RNCS()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.RNCS()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)



    return desc.iloc[0,0]


def RPCS(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.RPCS()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.RPCS()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def TASA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.TASA()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.TASA()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def TPSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.TPSA()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.TPSA()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def RASA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)

    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.RASA()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.RASA()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def RPSA(mol, **kwargs):

    seed = 429647

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    if mol.GetNumConformers() < 1: # added in 3.1 to avoid errors
        smi = Chem.MolToSmiles(mol, isomericSmiles = False)
        mol = Chem.MolFromSmiles(smi)
        AllChem.EmbedMolecule(mol, randomSeed=seed)
        if mol.GetNumConformers() < 1:
            desc = pd.DataFrame([np.nan])
        else:
            calc = Calculator(
                CPSA.RPSA()
                )
            desc = calc.pandas([mol], nproc=1, quiet=True)
    else:
        calc = Calculator(
            CPSA.RPSA()
            )
        desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


######### EState descriptors #########

def EStates(mol, **kwargs):
    calc = Calculator(EState.AtomTypeEState(kwargs['type'], kwargs['estate']))
    desc = calc.pandas([mol], nproc=1, quiet=True)

    return desc.iloc[0,0]


######### SLogP descriptors #########

def SLogP_PA(mol, **kwargs):


    calc = Calculator(
        SLogP.SLogP()
        )

    desc = calc.pandas([mol], nproc=1, quiet=True)


    return desc.iloc[0,0]


def SMR(mol, **kwargs):

    calc = Calculator(
        SLogP.SMR()
        )
    desc = calc.pandas([mol], nproc=1, quiet=True)

    return desc.iloc[0,0]



##############################################################################

def MoeType_all(mol):

    moetype_descriptors = {
        'LabuteASA': LabuteASA,
        'PEOE_VSA': PEOE_VSA,
        'SMR_VSA': SMR_VSA,
        'SLogP_VSA': SLogP_VSA,
        'EState_VSA': EState_VSA,
        'VSA_EState': VSA_EState,
        }

    unique = ['LabuteASA']

    k_values = [i for i in range(0, 13)]

    desc_df = pd.DataFrame()

    for desc in moetype_descriptors.keys():

         if desc in unique:
            try:
                desc_name = '{}'.format(desc)
                desc_value = moetype_descriptors[desc](mol)
                desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

                desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

            except: continue

         else:
            for k in k_values:
                kwargs = {'k': k}
                try:
                        desc_name = '{}{}'.format(desc, k)
                        desc_value = moetype_descriptors[desc](mol, **kwargs)
                        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

                        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

                except: continue

    return desc_df


def CPSA_all(mol):

    cpsa_descriptors = {
        'PNSA': PNSA,
        'PPSA': PPSA,
        'DPSA': DPSA,
        'FNSA': FNSA,
        'FPSA': FPSA,
        'WNSA': WNSA,
        'WPSA': WPSA,
        'RNCG': RNCG,
        'RPCG': RPCG,
        'RNCS': RNCS,
        'RPCS': RPCS,
        'TASA': TASA,
        'TPSA': TPSA,
        'RASA': RASA,
        'RPSA': RPSA,
        }

    uniques = ['RNCG', 'RPCG', 'RNCS', 'RPCS', 'TASA', 'TPSA', 'RASA', 'RPSA']

    version_values = [i for i in range(1, 6)]

    desc_df = pd.DataFrame()
    

    for desc in cpsa_descriptors.keys():

        if desc in uniques:
            try:
                desc_name = '{}'.format(desc)
                desc_value = cpsa_descriptors[desc](mol)
                desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

                desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
            except: continue

        else:
            for version in version_values:
                kwargs = {'version': version}
                try:
                    desc_name = '{}{}'.format(desc, version)
                    desc_value = cpsa_descriptors[desc](mol, **kwargs)
                    desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

                    desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
                except: continue

    return desc_df


def EState_all(mol):
    '''
    types:
        count --> N
        Sum --> S
        MAX --> MAX
        MIN --> MIN
    '''
    aggr_types = {
        'N' : 'count' ,
        'S' : 'sum',
        'MAX' : 'max',
        'MIN' : 'min'
        }

    es_types = ['sLi', 'ssBe', 'ssssBe', 'ssBH', 'sssB', 'ssssB', 'sCH3', 'dCH2',
    'ssCH2', 'tCH', 'dsCH', 'aaCH', 'sssCH', 'ddC', 'tsC', 'dssC', 'aasC', 'aaaC',
    'ssssC', 'sNH3', 'sNH2', 'ssNH2', 'dNH', 'ssNH', 'aaNH', 'tN', 'sssNH', 'dsN',
    'aaN', 'sssN', 'ddsN', 'aasN', 'ssssN', 'sOH', 'dO', 'ssO', 'aaO', 'sF', 'sSiH3'
    , 'ssSiH2', 'sssSiH', 'ssssSi', 'sPH2', 'ssPH', 'sssP', 'dsssP', 'sssssP', 'sSH'
    , 'dS', 'ssS', 'aaS', 'dssS', 'ddssS', 'sCl', 'sGeH3', 'ssGeH2', 'sssGeH',
    'ssssGe', 'sAsH2', 'ssAsH', 'sssAs', 'sssdAs', 'sssssAs', 'sSeH', 'dSe', 'ssSe'
    , 'aaSe', 'dssSe', 'ddssSe', 'sBr', 'sSnH3', 'ssSnH2', 'sssSnH', 'ssssSn', 'sI'
    , 'sPbH3', 'ssPbH2', 'sssPbH', 'ssssPb']

    desc_df = pd.DataFrame()

    for type, key in aggr_types.items():
        for estate in es_types:

            try:
                kwargs = {'type': key, 'estate': estate}
                desc_name = '{}{}'.format(type, estate)
                desc_value = EStates(mol, **kwargs)
                desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

                desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

            except: continue

    return desc_df

def SLogP_all(mol):

    slogp_descriptors = {
        'SLogP': SLogP_PA,
        'SMR': SMR,
        }

    desc_df = pd.DataFrame()

    for desc in slogp_descriptors.keys():

        try:
            desc_name = '{}'.format(desc)
            desc_value = slogp_descriptors[desc](mol)
            desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})
            desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)
        except: continue


    return desc_df
