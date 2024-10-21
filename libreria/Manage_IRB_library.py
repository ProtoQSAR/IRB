# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 10:43:10 2024

@author: proto (Eva Serrano-Candelas)

This script is used to manage the IRB library. It reads the sdf files and obtain the SMILES, setting a unique ID for each compound and 
merging all the files in a single one. It also creates a file with the unique SMILES and another with the duplicated SMILES.

The modifications are done in the following steps: Add the name variable = 'IRB_library_' to the output files, and generate the files for HYGIEIA, adding the column of unique "y".

modified by: Laureano E. Carpio (Moldrug)
"""

############################### CLEAN EVERYTHING #############################
for i in list(globals().keys()):
    if not i.startswith('_'):
        exec('del ' + i)
##############################################################################
#################################### IMPORTS #################################
import os
from pathlib import Path
import codecs
from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
from subprocess import call
import numpy as np
import pandas as pd
''' requires to install pip install xlrd '''
##############################################################################
################################ INITIAL VARIABLES ###########################
parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

data_folder =  '.' + os.path.sep + 'files'

results_folder = '.' + os.path.sep + 'pre_processing'

if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("Output directory created: ../results/")
else:
    print("Output directory already exists: ../results/")

name = 'IRB_library_'
###############################  FUNCTIONS ###################################

def retrieve_smiles(mol):

    return Chem.MolToSmiles(mol)

##############################################################################

############################## file 1 processing #############################

file1 = '2021 LIB_47489 CMPDS_26092024.sdf'

print(f'\n[+] Opening "{file1}" file')

df1 = PandasTools.LoadSDF(data_folder + os.path.sep + file1)

df1['SMILES']  = df1.apply(lambda row: retrieve_smiles(row['ROMol']), axis=1)

df1.insert(0, 'ID_SET', ['set_1']*df1.shape[0])

df1.index.name = 'ID_SET_COMPOUND'

output_file1_name = file1.split('.')[0]
df1.to_csv(results_folder+ os.path.sep + name + output_file1_name + '.csv', sep = ';')

##############################################################################

############################## file 2 processing #############################


file2 = 'all new library_104017 CMPDS_26092024.sdf'

print(f'\n[+] Opening "{file2}" file')

df2 = PandasTools.LoadSDF(data_folder + os.path.sep + file2)

df2.rename(columns = {'IRB PLATE ': 'IRB PLATE' }, inplace = True)

df2['SMILES']  = df2.apply(lambda row: retrieve_smiles(row['ROMol']), axis=1)

df2.insert(0, 'ID_SET', ['set_2']*df2.shape[0])

df2.index.name = 'ID_SET_COMPOUND'


output_file2_name = file2.split('.')[0]
df2.to_csv(results_folder+ os.path.sep + name + output_file2_name + '.csv', sep = ';')

##############################################################################

############################## file 3 processing #############################


file3 = 'focus 2 mM_70 cmpds_29092024_sent.sdf'

print(f'\n[+] Opening "{file3}" file')

df3 = PandasTools.LoadSDF(data_folder + os.path.sep + file3)

df3.rename(columns = {'IRB plate': 'IRB PLATE', 'Library': 'LIBRARY'}, inplace = True)

df3['SMILES']  = df3.apply(lambda row: retrieve_smiles(row['ROMol']), axis=1)

df3.insert(0, 'ID_SET', ['set_3']*df3.shape[0])

df3.index.name = 'ID_SET_COMPOUND'


output_file3_name = file3.split('.')[0]
df3.to_csv(results_folder+ os.path.sep + name + output_file3_name + '.csv', sep = ';')

##############################################################################

############################## file 4 processing #############################


file4 = 'focus 10 mM_10665 cmpds_29092024_sent.sdf'

print(f'\n[+] Opening "{file4}" file')

df4 = PandasTools.LoadSDF(data_folder + os.path.sep + file4)

df4.rename(columns = {'IRB plate ': 'IRB PLATE', 'Library': 'LIBRARY'}, inplace = True)

df4['SMILES']  = df4.apply(lambda row: retrieve_smiles(row['ROMol']), axis=1)

df4.insert(0, 'ID_SET', ['set_4']*df4.shape[0])

df4.index.name = 'ID_SET_COMPOUND'


output_file4_name = file4.split('.')[0]
df4.to_csv(results_folder+ os.path.sep + name + output_file4_name + '.csv', sep = ';')

##############################################################################

#%%

################################### merge ####################################


merged = pd.concat([df1,df2,df3,df4], axis = 0)

merged.reset_index(drop = True, inplace = True)

merged['ID_UNIQUE']  = 'IRB_' + merged.index.astype(str) 

merged.set_index('ID_UNIQUE', inplace = True)

merged_reordered = merged[['ID_SET', 'IRB WELL', 'LIBRARY', 'IRB PLATE', 'ID NUMBER', 'ID', 'ORIGINAL WELL', 'ROMol', 'SMILES']]

merged_reordered.to_csv(results_folder+ os.path.sep + name + 'merged_all.csv', sep = ';')

##################################prepare for HYGIEIA##########################

merged_hygieia = merged_reordered.copy()
merged_hygieia['y'] = 1
merged_hygieia.to_csv(results_folder+ os.path.sep + name + 'merged_all-preprocessed.csv', sep = ';')


##############################################################################

merged_unique_smiles = merged_reordered[~merged_reordered.duplicated(subset=['SMILES'], keep = False)]

merged_unique_smiles.to_csv(results_folder+ os.path.sep + name+ 'merged_all_unique_smiles.csv', sep = ';')

merged_duplicated_smiles = merged_reordered[merged_reordered.duplicated(subset=['SMILES'], keep = False)]

merged_duplicated_smiles.sort_values(by=['SMILES'], inplace = True)

merged_duplicated_smiles.to_csv(results_folder+ os.path.sep + name + 'merged_all_duplicated_smiles.csv', sep = ';')

print(f'\n[+] Files created in "{results_folder}" folder')

#print the number of molecules in the different datasets created
print(f'\n[+] Number of molecules in the different datasets created:')
print(f'File 1: {df1.shape[0]}')
print(f'File 2: {df2.shape[0]}')
print(f'File 3: {df3.shape[0]}')
print(f'File 4: {df4.shape[0]}')
print(f'Merged file: {merged_reordered.shape[0]}')
print(f'Merged file for HYGIEIA: {merged_hygieia.shape[0]}')
print(f'Unique SMILES: {merged_unique_smiles.shape[0]}')
print(f'Duplicated SMILES: {merged_duplicated_smiles.shape[0]}')

##################

##############END OF THE SCRIPT###############################################







