# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28

@author: Laureano E. Carpio (Moldrug)

This script is used to perform different checkpoints in the library sdfs. It reads the SDF and then check if the column "Mol" has "Empty mol" 
and if this is in line with the column "Cluster" that should be None.
"""
# Imports
from rdkit import Chem
import pandas as pd
from rdkit.Chem import PandasTools
import os
from pathlib import Path


################################ FUNCTIONS ###########################


def validate_empty_mol_condition(dataframe):
    """Verifica si se cumple la condición: 'Mol' == 'empty mol' implica 'Cluster' == None"""
    invalid_rows = dataframe[(dataframe["flag"] == "empty mol") & (dataframe["Cluster"] != "-")]
    if invalid_rows.empty:
        print("     [++] La condición de que si en flag aparece emtpy mol en Cluster aparece - se cumple en todos los casos.")
    else:
        print("     [++] La condición de que si en flag aparece emtpy mol en Cluster aparece - NO se cumple en las siguientes filas:")
        print(invalid_rows)


################################ INITIAL VARIABLES ###########################
parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

data_folder =  '.' + os.path.sep + 'clusterized'

############################## file 1 processing #############################

file1 = '2021 LIB_47489 CMPDS_26092024_clustered_clustered.sdf'

print(f'\n[+] Opening "{file1}" file')

df1 = PandasTools.LoadSDF(data_folder + os.path.sep + file1)

# checkpint
validate_empty_mol_condition(df1)

############################## file 2 processing #############################

file2 = 'all new library_104017 CMPDS_26092024_clustered_clustered.sdf'

print(f'\n[+] Opening "{file2}" file')

df2 = PandasTools.LoadSDF(data_folder + os.path.sep + file2)

# checkpint
validate_empty_mol_condition(df2)

############################## file 3 processing #############################

file3 = 'focus 2 mM_70 cmpds_29092024_sent_clustered_clustered.sdf'

print(f'\n[+] Opening "{file3}" file')

df3 = PandasTools.LoadSDF(data_folder + os.path.sep + file3)

# checkpint
validate_empty_mol_condition(df3)

############################## file 4 processing #############################

file4 = 'focus 10 mM_10665 cmpds_29092024_sent_clustered_clustered.sdf'

print(f'\n[+] Opening "{file4}" file')

df4 = PandasTools.LoadSDF(data_folder + os.path.sep + file4)

# checkpint
validate_empty_mol_condition(df4)