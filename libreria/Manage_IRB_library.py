# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 10:43:10 2024

@author: proto
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


###############################  FUNCTIONS ###################################

def retrieve_smiles(mol):
    
    return Chem.MolToSmiles(mol)


   

##############################################################################
    
############################## file 1 processing #############################
file1 = '2021 LIB_47489 CMPDS_26092024.sdf'

print(f'\n[+] Opening "{file1}" file')
    
df1 = PandasTools.LoadSDF(data_folder + os.path.sep + file1)

df1.replace('IRB PLATE ')

df1['SMILES']  = df1.apply(lambda row: retrieve_smiles(row['ROMol']), axis=1)

df1.insert(0, 'ID_SET', ['set_1']*df1.shape[0])

df1.index.name = 'ID_SET_COMPOUND'

output_file1_name = file1.split('.')[0]
df1.to_csv(results_folder+ os.path.sep + output_file1_name + '.csv', sep = ';')

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
df2.to_csv(results_folder+ os.path.sep + output_file2_name + '.csv', sep = ';')

##############################################################################
 
#%%
   
################################### merge ####################################


merged = pd.concat([df1,df2], axis = 0)

merged.reset_index(drop = True, inplace = True)

merged['ID_UNIQUE']  = 'IRB_' + merged.index.astype(str) 

merged.set_index('ID_UNIQUE', inplace = True)

merged_reordered = merged[['ID_SET', 'IRB WELL', 'LIBRARY', 'IRB PLATE', 'ID NUMBER', 'ID', 'ORIGINAL WELL', 'ROMol', 'SMILES']]

merged_reordered.to_csv(results_folder+ os.path.sep + 'merged_1_and_2.csv', sep = ';')



##############################################################################

merged_unique_smiles = merged_reordered[~merged_reordered.duplicated(subset=['SMILES'], keep = False)]



merged_unique_smiles.to_csv(results_folder+ os.path.sep + 'merged_1_and_2_unique_smiles.csv', sep = ';')

merged_duplicated_smiles = merged_reordered[merged_reordered.duplicated(subset=['SMILES'], keep = False)]

merged_duplicated_smiles.sort_values(by=['SMILES'], inplace = True)

merged_duplicated_smiles.to_csv(results_folder+ os.path.sep + 'merged_1_and_2_duplicated_smiles.csv', sep = ';')








