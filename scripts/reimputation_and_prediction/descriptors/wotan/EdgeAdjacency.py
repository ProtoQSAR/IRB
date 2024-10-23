import numpy as np
import pandas as pd

from rdkit import Chem

from math import log
from descriptors.utils.chem import EdgeAdjacencyMatrix
from utils.json_files import load_json


CONFIG = load_json('descriptors/config/descriptors.json')


def EdgeAdjacencyIndex(mol, **kwargs):
    '''Molecular descriptors calculated from the edge adjacency matrix of
    a molecule

    Keyword arguments:
        type (string) -- either:
            'conn': edge connectivity index
            'eigen': eigenvalues from edge adjancency matrix
            'spectral': spectral moment from edge adjacency matrix

        weight (string) -- property to weight on:
            'none': no weight
            'degrees': edge degrees
            'dipole': dipole moments
            'resonance': resonance integrals

        order (int) -- order of any of the three types
            [0, 1] for edge connectivity indices
            [1, 15] for eigenvalues and spectral moments

    'values for resonances from Dragon 5.5 manual'
    '''

    mol = Chem.RemoveHs(mol) # added in v3.0
    if kwargs['type'] == 'conn':
        raise NotImplementedError

    edge_adj_matrix = EdgeAdjacencyMatrix(mol)
        # changed in v3.0. It calls a external function

    bond_number = mol.GetNumBonds()
    bonds = mol.GetBonds()

    # Weight the edge adjacency matrix dependending on the selected property
    # NOTE: for 'degree', the sum of the rows of the edge adjacency matrix is
    # performed, which requires the matrix to be filled in the previous loop
    if kwargs['weight'] != 'none':
        # TODO: check correct SMARTS is retrieved from CONFIG
        # if kwargs['weight'] != 'degree':



        if kwargs['weight'] == 'bondorder':
            ''' Not implemented, but works well. v3.0
            Nevertheless, rdkit consider NO2 (nitro group) to have one single and
            one double bond. It gives different values thatn dragon7, guess that
            they consider two aromatic (1.5) bonds
            '''
            bond_orders = {
                    Chem.BondType.SINGLE: 1,
                    Chem.BondType.DOUBLE: 2,
                    Chem.BondType.TRIPLE: 3,
                    Chem.BondType.AROMATIC: 1.5
                }

            for j in range(len(edge_adj_matrix)):
                col = edge_adj_matrix[:,j]
                bond_order = bond_orders[bonds[j].GetBondType()]
                col[col != 0] = bond_order

        elif kwargs['weight'] == 'degree':
            ''' fixed in v3.0. It was computed erroneously and it didn't used
            weigthed matrices
            # edge_adj_matrix[i][i] = np.sum(edge_adj_matrix[:, i])
            '''
            for j in range(len(edge_adj_matrix)):
                col = edge_adj_matrix[:,j]
                edge_degree = np.sum(col)
                col[col != 0] = edge_degree

        elif kwargs['weight'] == 'resonance':
            resonance = CONFIG[kwargs['weight']]
            for j in range(len(edge_adj_matrix)):
                col = edge_adj_matrix[:,j]
                bond = bonds[j]

                bond_type = bond.GetBondType()

                atom1 = bond.GetBeginAtomIdx()
                atom1_symbol = mol.GetAtomWithIdx(atom1).GetAtomicNum()

                atom2 = bond.GetEndAtomIdx()
                atom2_symbol = mol.GetAtomWithIdx(atom2).GetAtomicNum()
                if [atom1_symbol,atom2_symbol] == [6, 8] or [atom1_symbol,atom2_symbol] == [8, 6]:
                    if str(bond_type) == 'DOUBLE':
                        resonance_value = 1.2
                    else:
                        resonance_value = 0.8
                elif str([atom1_symbol,atom2_symbol]) in resonance.keys():
                    key = str([atom1_symbol,atom2_symbol])
                    resonance_value = resonance[key]
                else:
                    resonance_value = 1
                col[col != 0] = resonance_value

        elif kwargs['weight'] == 'dipole':
            dipole = CONFIG[kwargs['weight']]
            idx_dict = dict()
            for key in CONFIG[kwargs['weight']].keys():
                smarts = Chem.MolFromSmarts(key)
                idx_dict[key] = sorted(mol.GetSubstructMatches(smarts))


            for key in CONFIG[kwargs['weight']].keys():
                smarts = Chem.MolFromSmarts(key)
                idx_dict[key] = sorted(mol.GetSubstructMatches(smarts))
            for j in range(len(edge_adj_matrix)):
                col = edge_adj_matrix[:,j]
                a1, a2 = (
                    mol.GetBondWithIdx(j).GetBeginAtomIdx(),
                    mol.GetBondWithIdx(j).GetEndAtomIdx()
                )

                for key in idx_dict.keys():

                    if (a1, a2) in idx_dict[key] or (a2, a1) in idx_dict[key]:
                        # Assume both atoms are the first and last atoms in key
                        col[col != 0] = CONFIG[kwargs['weight']][key]
            #     for j in range(len(edge_adj_matrix)):
            #         col = edge_adj_matrix[:,j]
            #         edge_degree = np.sum(col)
            #         col[col != 0] = edge_degree
# elif tuple([atom1_symbol,atom2_symbol]) in resonance.keys():
#     resonance_value = resonance[atom1_symbol,atom2_symbol]

    if kwargs['type'] == 'spectral':
        matrix_order = np.linalg.matrix_power(edge_adj_matrix, kwargs['order'])
        trace = np.trace(matrix_order)
        trace = float(trace)
        result = np.log(1 + trace)


    if kwargs['type'] == 'eigen':
        if not edge_adj_matrix.any(): return 0

        eigenvals = np.linalg.eigvalsh(edge_adj_matrix.tolist())
        if len(eigenvals) < kwargs['order']: return 0

        result =  eigenvals[-(kwargs['order'])]

    return result


def EdgeAdjacency_all(mol):

    desc_df = pd.DataFrame()

    dragon_name = {'spectral': 'ESpm', 'eigen': 'EEig'}
    weight_keys = {'none': 'u', 'degree': 'x', 'dipole': 'd', 'resonance': 'r'}
    # for type in ['spectral']:
    #     for weight in ['none', 'degree', 'dipole', 'resonance']:
    for type in ['spectral', 'eigen']:
        for weight in ['none', 'degree', 'dipole', 'resonance']:
            for order in range(1, 16):
                kwargs = {'type': type, 'order': order, 'weight': weight}
                desc_name = '{}{:02d}{}'.format(
                    dragon_name[type], order, weight_keys[weight]
                )
                desc_value = EdgeAdjacencyIndex(mol, **kwargs)
                desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

                desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

    return desc_df
