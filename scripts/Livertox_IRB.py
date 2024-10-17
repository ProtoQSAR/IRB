# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 13:29:55 2024

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
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
''' requires to install pip install xlrd '''
##############################################################################
################################ INITIAL VARIABLES ###########################
parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

input_path = '../datasets_for_modelling'

for (root, dirs, files) in os.walk(input_path):
    continue
#%%
# dict_properties = {}

# for file in files:
#     prop = file.split('_')[0]
#     if prop not in unique_properties:
#         unique_properties.append(prop)

# unique_properties = []

#%%

sdfs_to_process = ['Oatp1b1INH_test.sdf','Oatp1b1INH_training.sdf','Oatp1b3INH_test.sdf','Oatp1b3INH_training.sdf']


dict_properties = {}

print('[+] Processing files')

for sdf in files:
    
   
    if sdf in sdfs_to_process:
        print("\t", sdf)
        
        prop = sdf.split('_')[0]
        
        if prop not in dict_properties.keys():
            dict_properties[prop] = {}
        
        dataset = PandasTools.LoadSDF(input_path + os.path.sep + sdf)
        
        smiles = [Chem.MolToSmiles(mol) for mol in dataset['ROMol']]
        
        dataset.rename(columns={'SMILES': 'orig_SMILES'}, inplace = True)
        
        dataset.insert(0,'SMILES',smiles)
        
        sdf.split('.')[0] = sdf.split('_')[1]
        
        subset = sdf.split('.')[0].split('_')[1]
        
        list_subset = dataset.shape[0]*[subset]
        
        dataset.insert(dataset.shape[1],'subset',list_subset)
        
        dataset = dataset.drop(['ROMol'], axis='columns')
        
        dataset.index.name = 'idx_separated'
        
        dict_properties[prop][subset] = dataset
        
        dict_properties[prop][subset + '_size'] = dataset.shape
        
        dict_properties[prop][subset + '_columns'] = list(dataset.columns)
        outputfilename = sdf.split('.')[0]
        
        # dataset.to_csv(f'./pre_preprocessed_ONGOING/input_data/Livertox_{outputfilename}.csv', index = True, sep = ';')

#%%
print('\n[+] Merging files')

for key, value in dict_properties.items():
    
    
    
    if 'test' in value.keys():
        df_training = value['training']
        df_test = value['test']

        df_test.rename(columns = {'Binary Characterization':'Binary_Characterization',
                                  'Binary characterization':'Binary_Characterization'}, inplace = True)
        
        
        df_all = pd.concat([df_training, df_test], axis = 0)
        dict_properties[key]['all'] = df_all
        dict_properties[key]['all_size'] = df_all.shape
    
        # print('\t\ttrain: ', dict_properties[key]['training_size'], ', ', df_training.shape)
        # print('\t\ttest: ', dict_properties[key]['test_size'], ', ',  df_test.shape)
        # print('\t\tall_prev: ', dict_properties[key]['training_size'][0] + dict_properties[key]['test_size'][0])
        # print('\t\tall_now: ', df_training.shape[0] + df_test.shape[0], ', ',  df_all.shape[0])
        
        df_all.reset_index(inplace = True)
        
        df_all.index.name = 'idx_together'
        
        df_all.to_csv(f'./pre_preprocessed_ONGOING/input_data/Livertox_{key}_all.csv', index = True, sep = ';')
    
    else:
        df_all = value['training']
        # print('\t\ttrain: ', dict_properties[key]['training_size'], ', ', df_training.shape)
        
        df_all.to_csv(f'./pre_preprocessed_ONGOING/input_data/Livertox_{key}_only_training.csv', index = True, sep = ';')
        
    if 'TRANSP' in key:
        print("\t", key)
        print('\t\t negative', df_all['Activity'].value_counts()[0])
        print('\t\t positive', df_all['Activity'].value_counts()[1])
        
        
        
