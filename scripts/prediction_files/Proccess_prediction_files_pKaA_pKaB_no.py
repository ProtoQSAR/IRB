# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 09:41:29 2024

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
import pandas as pd
import numpy as np
from rdkit import Chem
##############################################################################
################################ INITIAL VARIABLES ###########################

pd.set_option('mode.chained_assignment', None) # pd.set_option('mode.chained_assignment', 'warn')

parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')
data_folder =  '..' + os.path.sep + 'compare_ONGOING' + os.path.sep + 'results'


##############################################################################

###############################  FUNCTIONS ###################################



def create_dir_for_io(folder):

    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"\t\tDirectory created: {folder}")
    else:
        print(f"\t\tDirectory already exists: {folder}")
        
def sanitize_smiles(row):
    mol = Chem.MolFromSmiles(row)
    
   
    san_smi = Chem.MolToSmiles(mol, isomericSmiles = False)

    return san_smi
        
##############################################################################

############################### CONFIGURATION FILES ##########################        

data_folder =  '..' + os.path.sep + 'compare_ONGOING' + os.path.sep + 'results'

root_path = '.' + os.path.sep + 'prediction_files'  + os.path.sep + 'pKa' 
        
prediction_folder = root_path + os.path.sep + 'predictions'

originals_folder = root_path + os.path.sep + 'originals'

merging_folder = root_path + os.path.sep + 'merged'

create_dir_for_io(merging_folder)



##############################################################################




endpoints = ['pKaA', 'pKaB']

endpoint= 'pKaA'

datasets = ['avdeef', 'sander', 'liao', 'settimo']

datasets = ['avdeef', 'sander']



if endpoint == 'pKaA':
    tag1 = 'acidic'
    tag2 = 'Acid'

elif endpoint == 'pKaB':
    tag1 = 'basic'
    tag2 = 'Basic'


df_all0 = pd.read_csv(data_folder + os.path.sep + f'all_merged_PhCh_{endpoint}.csv', sep = ';')
# df_all_reindexed = df_all0.set_index('ID_all')

for dataset in datasets:
    
    # df_original = pd.read_excel(originals_folder + os.path.sep + f'PhCh_pKa_{dataset}_for{tag2}.xlsx', sheet_name= f'pka_{dataset}')
    
    df_all0['forindexinall'] = df_all0[f'orig_ID_{dataset}{tag2}']
    
    df_all0.set_index('forindexinall', inplace = True)
    
    
    df_predicted = pd.read_excel(prediction_folder + os.path.sep + f'pka_{tag1}.xlsx', sheet_name= f'pka_{dataset}')
    
   

    
    
    df_predicted.rename(columns = {'ID': f'ID_{dataset}',
                                   'SMILES': f'SMILES_{dataset}',
                                   'Type': f'Type_{dataset}',
                                   'type': f'Type_{dataset}',
                                   'pKa_Exp': f'pKa_Exp_{dataset}',
                                   
                                   'OPERA_pKa_a_Pred': 'OPERA_pKa-a_Pred',
                                   'OPERA_pKa_b_Pred': 'OPERA_pKa-b_Pred',
                                   

                                   }, inplace = True)
    

    
        
    df_predicted[f'SMILES_{dataset}noniso'] = df_predicted.apply(lambda row: sanitize_smiles(row[f'SMILES_{dataset}']), axis=1)

    df_predicted = df_predicted[['ChemAxon_pKa-a_Pred', 
                                'ChemAxon_pKa-b_Pred', 'ChemAxon_pKa_Pred',
                                'OPERA_pKa-a_Pred', 'OPERA_pKa-b_Pred', 
                                'OPERA_pKa_Pred', 'OPERA_pKa_AD', 
                                'OPERA_pKa_TS', f'ID_{dataset}', f'SMILES_{dataset}', f'Type_{dataset}', 
                                f'pKa_Exp_{dataset}', f'SMILES_{dataset}noniso']    ]






    df_predicted['forindexinpred'] = df_predicted[f'ID_{dataset}']
    
    df_predicted.set_index('forindexinpred', inplace = True)
    
    df_all0['SMILES_noniso'] = df_all0.apply(lambda row: sanitize_smiles(row['SMILES']), axis=1)
    
    
    
    diff = df_predicted[~df_predicted[f'SMILES_{dataset}noniso'].isin(df_all0['SMILES_noniso'])]
    
    if diff.shape[0] == 0:
        print('all mols found :) ')
    else:
        print('not a lucky girl, sorry :(')
        
    print(df_predicted.shape)
    print(list(df_predicted.columns))
    
    
    
    print(df_predicted.head())
    
    print(df_all0.columns)
    
    df_all0 = pd.concat([df_all0,df_predicted], axis = 1, join="outer")
    
    
    
    
    
    
    # df_all.set_index()
    
#     smiles.append(list(df_predicted['SMILES_noniso']))
        
        
#     df_predicted.set_index('SMILES_noniso', inplace = True)
        
#     list_df_predicted.append(df_predicted)
    
#     num_mols = num_mols + df_predicted.shape[0]
    
    
    
# print(num_mols)



# flattenlist = [x for xs in smiles for x in xs]

# smiles_as_set = set(flattenlist)

# print(len(smiles_as_set))

# #%%
# df_0 = list_df_predicted[0]

# for df in list_df_predicted[1:]:

#     # df_0 = pd.concat([df_0,df], axis = 1)
    
#     df_0 = df_0.append(df)
    
# df_0.to_csv(merging_folder + os.path.sep + 'mii.csv', sep = ';', index = True)
    
# # merged = pd.concat(list_df_predicted,axis = 0)

# df_final = df_0.sort_values('SMILES_noniso').groupby('SMILES_noniso').apply(lambda x: x.ffill())

# df_final = df_final.sort_values('SMILES_noniso').groupby('SMILES_noniso').apply(lambda x: x.bfill())

# df_final.to_csv(merging_folder + os.path.sep + 'miooo.csv', sep = ';', index = True)


# #%%
# smi = 'CC(=CCl)C(=O)O'

# mol = Chem.MolFromSmiles(smi)

# san_smi = Chem.MolToSmiles(mol, isomericSmiles = False)

# print(san_smi)



















