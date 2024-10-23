
from rdkit import Chem
import os
from itertools import chain

from utils.json_files import load_json

import pandas as pd

all_dicts = load_json('descriptors/config/Atom_centered_v3_0.json')


def atom_centred(mol, **kwargs):
    # Get specified descriptor
    #print(kwargs) 
    smarts = [Chem.MolFromSmarts(sm)
              for sm in all_dicts["smarts"][kwargs['key']]]
    # Get position of central atom in the tuple generated  by GetSustructMatches
    posi_atom_centr = [
        posi for posi in all_dicts["central_atoms"][kwargs['key']]]

    # count the matchesof the fragment defined by the descriptor in the molecule
    total_matches = []
    # Iterate through the list of smarts and get "num"which is the number of the code in the list
    for num, sm in enumerate(smarts):
        # gets a  list of tuples that identify the matches in the molecule
        mol = Chem.AddHs(mol)
        matches = mol.GetSubstructMatches(sm, uniquify=True)
        # get the number that identifies the position of the central atom in the fragment
        posi = posi_atom_centr[num]
        # remove repetitions (for a certain central atom, only count one match)
        norep_matches = []
        for i, match in enumerate(matches):
            if i == 0:
                norep_matches.append(matches[i])
            elif matches[i][posi] != matches[i-1][posi]:
                norep_matches.append(matches[i])

        total_matches.append(len(norep_matches))

    return sum(total_matches)


def atom_centred_hybr(mol, **kwargs):
    # Get specified descriptor
    #print(kwargs) 
    smarts = [Chem.MolFromSmarts(sm)
              for sm in all_dicts["smarts"][kwargs['key']]]
    mol = Chem.AddHs(mol)
    # Get which is the central atom of the descriptor

    posi_atom_centr = [
        posi for posi in all_dicts["central_atoms"][kwargs['key']]]
    # Get which is the central carbon atom of the descriptor
    c_frag = [posi for posi in all_dicts["c_attached_to_central"][kwargs['key']]]
    # Get which is the oxidation number required for the descriptor
    oxy_num = [
        num for num in all_dicts["oxydation_number_descriptor"][kwargs['key']]]

    # Get a list of all atoms in the molecule
    atom_list = [atom for atom in mol.GetAtoms()]

    total_matches = []
    # Go through the list of smarts to match for each descriptor (normally there is only one, but there are some taht have up to 4 different options)
    for num, sm in enumerate(smarts):
        if sm in all_dicts["heteroatoms_FClBrI"]:
            norep_matches2 = mol.GetSubstructMatches(sm, uniquify=True)
        else:
            # find the smarts in the molecule
            matches = mol.GetSubstructMatches(sm, uniquify=True)
            # get the central atom
            posi = posi_atom_centr[num]
            # get the main C attached to the central atom of the descriptor
            c_frag_num = c_frag[num]
            # get the oxidation number defined in the descriptor
            oxy_num_num = oxy_num[num]

            # get matches of the fragment in the molecule and remove repetitions with the same central atom
            norep_matches = []
            for i, match in enumerate(matches):
                if i == 0:
                    norep_matches.append(matches[i])
                elif matches[i][posi] != matches[i-1][posi]:
                    norep_matches.append(matches[i])
            # get the oxidation number of the C attached to the central atom and consider only
            # the fragment that has the correct oxidation number defined in the descriptor
            norep_matches2 = []
            for match in norep_matches:
                # isolate the carbon atom
                carbon = match[c_frag_num]
                carbonObj = atom_list[carbon]
                bonds = carbonObj.GetBonds()
                c_oxi = []
                # work with the bonds the C forms with adjacent atoms
                n_pir = False
                for bond in bonds:
              
                    c_oxi_ind = 0
                    
                    # if C#N, it is+2and not +3
                    #if str(bond.GetEndAtom().GetSymbol()) == "N" and str(bond.GetBondType()) == "TRIPLE":
                    if str(bond.GetOtherAtom(carbonObj).GetSymbol()) == "N" and str(bond.GetBondType()) == "TRIPLE":
                        c_oxi_ind = 2
                    else:
                        neighbor_value = all_dicts["negative_energies"][str(
                            #bond.GetEndAtom().GetSymbol())]
                            bond.GetOtherAtom(carbonObj).GetSymbol())]
                        bond_value = all_dicts["bond_types"][str(
                            bond.GetBondType())]
                        
                        c_oxi_ind = neighbor_value * bond_value
                        if str(bond.GetOtherAtom(carbonObj).GetSymbol()) == "N" and bond.GetOtherAtom(carbonObj).GetTotalDegree() == 2 and str(bond.GetBondType()) == "AROMATIC":
                            n_pir = True
                    c_oxi.append(c_oxi_ind)
                if n_pir == True:
                    c_oxi.append(1)
                oxi_total = sum(c_oxi)
                
                # oxidation number of 5 is the code given to identify the descriptor that must have
                # an oxidation metween 2 and 4
                if oxy_num_num == 5:
                    if oxi_total >= 2 and oxi_total <= 4:
                        norep_matches2.append(match)
                # there are 5 descriptors with H attached to C with oxidation number 0
                # that specify the number of heteroatoms attached to adjacent C atoms:
                elif oxy_num_num == 0 and kwargs['key'] in list(all_dicts["heteroatoms_descriptorsH"].keys()):
                    if oxi_total == 0:
                        count_x = 0
                        Neighbors = carbonObj.GetNeighbors()
                        for neighbor in Neighbors:
                            neighbor_symbol = str(neighbor.GetSymbol())
                            if neighbor_symbol == "C":
                                neineis = neighbor.GetNeighbors()
                                for neinei in neineis:
                                    if str(neinei.GetSymbol()) in all_dicts["heteroatoms"]:
                                        count_x = count_x+1
                        #print(count_x,all_dicts["heteroatoms_descriptorsH"][kwargs['key']])
                        if count_x >= 4 and all_dicts["heteroatoms_descriptorsH"][kwargs['key']] == 4:
                            norep_matches2.append(match)
                        elif count_x == all_dicts["heteroatoms_descriptorsH"][kwargs['key']]:
                            
                            norep_matches2.append(match)
                # if oxidaton number obtained for the C in the fragment is the oxidation number
                # required in the descriptor
                elif oxi_total == oxy_num_num:
                    #print(oxi_total == oxy_num_num,oxi_total,oxy_num_num)
                    norep_matches2.append(match)
        # count the real number of matches after clearing all repetitions
        total_matches.append(len(norep_matches2))
        #print(total_matches)
    return sum(total_matches)


def AtomCentred_all(mol):

    desc_df = pd.DataFrame()

    # transform SMILES to mol object
    mol = Chem.AddHs(mol)

    for key in all_dicts["smarts"]:

        kwargs = {'key': key}
        desc_name = key

        if key in all_dicts["descr_hybridation"]:
            desc_value = atom_centred_hybr(mol, **kwargs)
        else:
            desc_value = atom_centred(mol, **kwargs)

        desc_value = pd.DataFrame.from_dict({desc_name: [desc_value]})

        desc_df = pd.concat([desc_df, desc_value], axis=1, sort=False)

    return desc_df
