import re
import numpy as np

from rdkit import Chem
from rdkit.Chem import (
    rdmolops, Atom, PeriodicTable, Graphs, rdchem, rdMolDescriptors
)

from mordred._atomic_property import *

from networkx import Graph, floyd_warshall_numpy

from math import pi

from utils.json_files import load_json

CONFIG = load_json('descriptors/config/descriptors.json')


def load_descriptors(path):
    '''Loads and substitutes extrange characters in python in order to
    call functions
    '''
    with open(path, 'r') as infile:
        content = infile.read().split('\n')
        content = [re.sub('[\[\]-]', '', line) for line in content]

    return list(filter(None, content))


def replace_atoms_two_letters(smiles, letter='Q'):
    for atom in CONFIG['TWO_LETTERS']:
        smiles = smiles.replace(atom, letter)

    return smiles


def get_atoms(smiles):
    smiles = replace_atoms_two_letters(smiles)
    return re.sub(r'[^a-zA-Z]+', '', smiles)


def count_cycles(smiles, length):
    '''Counts the numbers of cycles and their length
    '''
    numbers = set([int(s) for s in smiles if s.isdigit()])

    # Get the cycles between numbers
    cycles = [re.search('%s(.*)%s' % (n, n) , smiles).group(1)
              for n in numbers]

    # Fix errors while dealing with atoms represented by more than one letter
    cycles = [replace_atoms_two_letters(cyc) for cyc in cycles]

    # Add one in length for the first atom before the first number
    cycles_len = [len(cyc) + 1 for cyc in cycles]

    # Filter those lengths according to the descriptor's definition
    cycles_len = [cyc_len for cyc_len in cycles_len
                  if cyc_len == length]

    return len(cycles_len)


def GetPrincipleQuantumNumber(atNum):
    '''Get principal quantum number for atom number
    '''
    if atNum <= 2:
        return 1
    elif atNum <= 10:
        return 2
    elif atNum <= 18:
        return 3
    elif atNum <= 36:
        return 4
    elif atNum <= 54:
        return 5
    elif atNum <= 86:
        return 6
    else:
        return 7


def get_atom_radii(atom_list):
    radii = list()
    # revised in v3.0. Added S3 to radii_list
    hybridization = []
    for atom in atom_list:
        symbol = atom.GetSymbol()
        valence = atom.GetHybridization()
        hybridization.append(str(atom.GetHybridization()))
        try:
            if symbol in ['C','N','O','P']:
                radius_key = '{}_{}'.format(symbol, valence)
            else:
                radius_key = '{}'.format(symbol)
            # print(radius_key)
            radii.append(CONFIG['radii_v3'][radius_key])
        except:
            radii.append(0)
    # print(hybridization)


    return radii


def volume_vdw(vdw_rad):
    vdw_vol = (4/3)*pi*vdw_rad**3
    return vdw_vol


def laplacian_eigenvalues(mol, nSK, nBO):

    # Get the vertex degree matrix
    bond_list = []
    for bond in range(nBO):
        bond_list.append(mol.GetBondWithIdx(bond).GetBeginAtomIdx())
        bond_list.append(mol.GetBondWithIdx(bond).GetEndAtomIdx())

    vertex_matrix = np.zeros(shape=(nSK, nSK), dtype=np.object)
    i = 0
    for i in range(nSK):
        vertex_matrix[i,i] = bond_list.count(i)

    # Get the adjacency matrix and calculate the laplace matrix
    adj_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol)

    laplace = np.subtract(vertex_matrix, adj_matrix)

    # Laplace matrix has to be converted into a array-list to get eigenvals
    eigenval = np.linalg.eigvalsh(laplace.tolist())

    return eigenval


def EdgeAdjacencyMatrix(mol):

    bond_number = mol.GetNumBonds()
    bond_list = []
    for bond in range(bond_number):
        bond_list.append((
            mol.GetBondWithIdx(bond).GetBeginAtomIdx(),
            mol.GetBondWithIdx(bond).GetEndAtomIdx()
        ))

    # WARNING: keep type of matrix as np.object. Otherwise, if the size of int
    # is reduced (such as np.int8), the exponentation will yield negative values
    # See: https://stackoverflow.com/questions/39602404/numpy-matrix-exponentiation-gives-negative-value
    edge_adj_matrix = np.zeros(shape=(bond_number, bond_number), dtype=np.object)

    # Fill the edge adjacency matrix
    for i in range(bond_number):
        a1, a2 = (
            mol.GetBondWithIdx(i).GetBeginAtomIdx(),
            mol.GetBondWithIdx(i).GetEndAtomIdx()
        )

        # Get the other pair of atoms in order to detect bond/edge adjancency
        for j in range(bond_number):
            a3, a4 = (
                mol.GetBondWithIdx(j).GetBeginAtomIdx(),
                mol.GetBondWithIdx(j).GetEndAtomIdx()
            )

            # If bonds next to each other, one atom is the same in both bonds
            if (a1 == a3 or a1 == a4 or a2 == a3 or a2 == a4) and i != j:
                edge_adj_matrix[i][j] = 1
                edge_adj_matrix[j][i] = 1
    return edge_adj_matrix


def detour_matrix(mol):
    '''The detour matrix (square symmetric matrix representing a H-depleted
       molecular graph, whose entry i-j is the length of the longest path from
       vertex vi to vertex vj) is used to compute three descriptors.
    '''
    mol = Chem.RemoveHs(mol)
    nSK = mol.GetNumAtoms()
    detour = np.zeros(shape=(nSK, nSK))

    # Calculate the detour matrix
    for distance in reversed(range(1, nSK)):
        for path in Chem.FindAllPathsOfLengthN(mol, distance+1, useBonds=0):
            # Sometimes it finds paths that cross an atom twice, so we want to
            # avoid this possibility checking if the last index is already in
            # the list of indexes of the path
            if [i for i in path].count(path[-1]) != 1:
                continue

            if detour[path[0], path[-1]] == 0 and path[0] != path[-1]:
                detour[path[0], path[-1]] = distance
                detour[path[-1], path[0]] = distance

    return detour


def GetBaryszDistanceMatrix(mol, **kwargs):
    ''' Get the Barysz distance matrix '''

















def GetWeightedMatrix(mol, **kwargs):
    '''The Wiener-type indices from weighted distance matrices (Whetw) are
       calculated by using the same formula as the Wiener index W applied to
       each weighted distance matrix, i.e. half-sum of matrix entries.

       Follows the mordred implementation of Barysz matrix # comment in v3.0.
       The results are different from Dragon
    '''
    prop_fn_dict = {
        'Z': lambda atom: atom.GetAtomicNum(),
        'm': lambda atom: atom.GetMass(),
        'v': get_vdw_volume,
        'e': get_sanderson_en,
        'p': get_polarizability,
    }
    prop_fn = prop_fn_dict[kwargs['weight']]


    # Helper functions
    get_weight = lambda atom_i, atom_j, pi, fn: \
        (C * C) / (fn(atom_i) * fn(atom_j) * pi)

    fill_diagonal = lambda sp, mol, fn: \
        np.fill_diagonal(sp, [1. - C / fn(a) for a in mol.GetAtoms()])


    # Initialize carbon property value
    C = prop_fn(Chem.Atom(6))

    G = Graph()

    G.add_nodes_from(a.GetIdx() for a in mol.GetAtoms())

    # it fills the nodes whete i !=j
    for bond in mol.GetBonds():

        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        atom_i = bond.GetBeginAtom()
        atom_j = bond.GetEndAtom()

        pi = bond.GetBondTypeAsDouble()

        w = get_weight(atom_i, atom_j, pi, prop_fn)

        G.add_edge(i, j, weight=w)


    sp = floyd_warshall_numpy(G)

    # it fills the nodes whete i ==j
    fill_diagonal(sp, mol, prop_fn)

    return sp


def BalabanJ(mol, dMat):
    # WARNING: adapted from RDkit
    '''Calculate Balaban's J value for a molecule

        **Arguments**
            - mol: a molecule
            - dMat: (optional) a distance/adjacency matrix for the molecule, if this
            is not provide, one will be calculated
            - forceDMat: (optional) if this is set, the distance/adjacency matrix
            will be recalculated regardless of whether or not _dMat_ is provided
            or the molecule already has one

        **Returns**
            - a float containing the J value

    We follow the notation of Balaban's paper:
    Chem. Phys. Lett. vol 89, 399-404, (1982)
    '''
    adjMat = Chem.GetAdjacencyMatrix(
        mol, useBO=0, emptyVal=0, force=0, prefix='NoBO'
    )

    # s = sum(dMat).tolist()[0] # changed in version v3.0. it returns a np.matrix instead of an nparray
    s = sum(dMat).tolist()

    q = mol.GetNumBonds()
    n = mol.GetNumAtoms()
    mu = q - n + 1

    suma = 0.
    nS = len(s)
    for i in range(nS):
        si = s[i]
        for j in range(i, nS):
            if adjMat[i, j] == 1:
                suma += 1. / float(np.sqrt(si * s[j]))

    if mu + 1 != 0:
        J = float(q) / float(mu + 1) * suma
    else:
        J = 0

    return J

#DEPRECATED
# give inconsistent results --> revised in v3.0 toValenceVertexDegree2
# def ValenceVertexDegree(atom):
#     "Computes the valence vertex degree of an atom"
#
#     qn = GetPrincipleQuantumNumber(atom.GetAtomicNum())
#
#     if atom.GetAtomicNum() == 7 and atom.GetIsAromatic():
#         Hs = rdchem.Atom.GetNumExplicitHs(atom)
#
#     else:
#         Hs = rdchem.Atom.GetNumImplicitHs(atom)
#
#     #Dictionary to relate the qn with the value to compute the valence vertex degree
#     qn_dict = {2:2, 3:10, 4:28, 5:46, 6:64}
#
#     return atom.GetAtomicNum() - qn_dict[qn] - Hs

#DEPRECATED
# definitive revision in version v3.0
def ValenceVertexDegree2(atom):
    "Computes the valence vertex degree of an atom"
    ''' Defined as in Maestro and old versions of Dragon (<5.5)
        (Z^v_i - h_i) / (Z_i - Z^v_i - 1) where:
        Z^v_i = number of valence electrons of ith atom
        h_i = number of hydrogens bonded to atom
        Z_i = total number of electrons of ith atom (atomic number)
    '''


    tb1 = Chem.GetPeriodicTable()
    Zv_i = tb1.GetNOuterElecs(atom.GetAtomicNum())
    # print(atom.GetAtomicNum(), Zv_i)
    h_i = atom.GetTotalNumHs()
    result = float(Zv_i - h_i)/(atom.GetAtomicNum() - Zv_i - 1)



    return result

# this is more similar to Dragon7, but gives problems in other descriptors
def ValenceVertexDegree3(atom):
    "Computes the valence vertex degree of an atom"
    ''' Defined as in upper versions of Dragon (>5.5) simplified formula
        vertex degree minus the number of attached atoms,
        formula taken from MATCH Commun. Math. Comput. Chem. 64(2010) 359-372
        Z^v_i - h_i
        Z^v_i = number of valence electrons of ith atom
        h_i = number of hydrogens bonded to atom
    Fixed by using Maestro documentation
    '''

    tb1 = Chem.GetPeriodicTable()
    Zv_i = tb1.GetNOuterElecs(atom.GetAtomicNum())

    h_i = atom.GetTotalNumHs()
    # print( '\n\tZv_i', Zv_i ,'\n\tatom.atomic_number', atom.GetAtomicNum(),'\n\th_i ', h_i)
    result = (Zv_i -atom.GetTotalNumHs())

    return result
