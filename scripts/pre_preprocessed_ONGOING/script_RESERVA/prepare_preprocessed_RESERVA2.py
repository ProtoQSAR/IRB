# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 18:45:54 2024

@author: proto
"""
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


##############################################################################
################################ INITIAL VARIABLES ###########################

parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

convert_ascii_data = True
data_folder =  '.' + os.path.sep + 'input_data'

results_folder = '.' + os.path.sep + 'results'

if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("Output directory created: ../results/")
else:
    print("Output directory already exists: ../results/")
    



path_to_higieia = r'F:\ProtoQSAR\GENERATE_MODELS_3_0\FullAsFunctions'
exec_hygieia = 'hygieia_mod.py'

config_file = 'info_original_datasets.xlsx'

##############################################################################

config_df = pd.read_excel('..' + os.path.sep + config_file)

config_df = config_df[config_df['proceed'] ==  'yes']

config_df_sel_cols = config_df[['dataset',
       'curation_process', 'sheet', 'encoding', 'experimental_column', 'units',
       'output_name', 'require_preprocessing', 'processed']]


dictios_dfs = config_df_sel_cols.set_index('dataset').T.to_dict()



''' 
Please, check the following dicts to be sure about renaming columns
'''

#%%


dictio_posibles_renames_smiles = {'smiles':'SMILES',
                                  'SMILES.desalt':'SMILES',
                                  'Smiles':'SMILES',
                                  'Smiles string':'SMILES',
                                  'SMILES structure': 'SMILES'}

dictio_posibles_renames_ids = {'ID' : 'orig_ID',
                               'Structure [idcode]' : 'orig_ID',
                               'Nº' : 'orig_ID',
                               '#' : 'orig_ID',
                               'id' : 'orig_ID',
                               'No': 'orig_ID'}

for dataset, info in dictios_dfs.items():

    print(f'\n[+] Analysing "{dataset}" dataset')
    
    extension = dataset.split('.')[1]
    
    pathfile = data_folder + os.path.sep + dataset
    
    file_convertedname = dataset.split('.')[0]
    pathfile_converted = data_folder + os.path.sep + f'{file_convertedname}-converted.sdf'
    
   
    if info['curation_process'] == 'curation1':   

        if info['encoding'] != 'utf8':
            
            print(f'\t[++] file is in {info["encoding"]} encoding and wil be converted to "utf8" encoding')
            
            #read input file
            with codecs.open(pathfile, 'r', encoding = info['encoding']) as file:
              lines = file.read()
            
            #write output file
            with codecs.open(pathfile_converted, 'w', encoding = 'utf8') as fileo:
              fileo.write(lines)
            
            
            print(f'\t\t[+++] Converted file created: {file_convertedname}_converted.sdf')
        
            df = PandasTools.LoadSDF(pathfile_converted)
        else:
            df = PandasTools.LoadSDF(pathfile)
            
        original_output_file_df = f'{info["output_name"]}-original.csv'
        print(f'\t[++] Dataframe file created: {original_output_file_df}')
        df.to_csv(results_folder + os.path.sep + original_output_file_df, sep = ';', index = False)  

                  
        if 'SMILES' in df.columns:
            df.rename(columns={'SMILES': 'orig_SMILES'}, inplace = True)
          
        
        smiles = [Chem.MolToSmiles(mol) for mol in df['ROMol']]
        df.insert(0,'SMILES',smiles)
            
            
    elif info['curation_process'] == 'curation2':
        
        if info['require_preprocessing'] == 'require_preprocessing_3':
        
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], header = 1)
        
        elif info['require_preprocessing'] == 'require_preprocessing_4':
            
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], header = 1)
        
        else:
            df = pd.read_excel(pathfile, sheet_name = info['sheet'])
        
        original_output_file_df = f'{info["output_name"]}-original.csv'
        print(f'\t[++] Dataframe file created: {original_output_file_df}')
        df.to_csv(results_folder + os.path.sep + original_output_file_df, sep = ';', index = False)
        
        
        if info['require_preprocessing'] == 'require_preprocessing_1':
            df.rename(columns={'SMILES': 'treated_SMILES'}, inplace = True)
            
        if info['require_preprocessing'] == 'require_preprocessing_2':
            allowed_endpoints = ['pKa', 'Acidic pKa', 'pKa(1)']
            df = df[df['Endpoint'].isin(allowed_endpoints)]
            df = df[df['Qualifier'] == 'No qualifier']
            
        if info['require_preprocessing'] == 'require_preprocessing_4':
            df = df[df['SDg'] < 2]
            

        
#%%        
        
    for column in list(df.columns):
        if column in dictio_posibles_renames_smiles.keys():
            
            df.rename(columns = dictio_posibles_renames_smiles, inplace = True)

    for column in list(df.columns):
        if column in dictio_posibles_renames_ids.keys():
            
            df.rename(columns = dictio_posibles_renames_ids, inplace = True) 

    if 'orig_ID' not in list(df.columns):
        df['orig_ID'] = [np.nan]*df.shape[0]
            
    
   
    df.index.name = 'ID'
    df.reset_index(inplace = True)
    
    df.insert(df.shape[1],'experimental_column',info["experimental_column"])
    df.insert(df.shape[1],'UNITS',info["units"])
    
    entire_output_file_df = f'{info["output_name"]}-entiredata.csv'
    print(f'\t[++] Dataframe file created: {entire_output_file_df}')
    df.to_csv(results_folder + os.path.sep + entire_output_file_df, sep = ';', index = False)
    
    
    
    if 'CAS' not in list(df.columns):
        df['CAS'] = [np.nan]*df.shape[0]
    if 'NAME' not in list(df.columns):
        df['NAME'] = [np.nan]*df.shape[0]
    else:
        df['NAME'] = df['NAME'].str.replace('|', '_', regex = True)
        df['NAME'] = df['NAME'].str.replace(';', ',', regex = True)  
        
    cols_to_maintain = ['ID', 'orig_ID', 'SMILES', info["experimental_column"], 'CAS', 'NAME', 'UNITS'] 
    
    df_preprocessed = df[cols_to_maintain]
    
    df_preprocessed.rename(columns={info["experimental_column"]: 'y'}, inplace = True)
   
    preprocessed_output_file_df = f'{info["output_name"]}-preprocessed.csv'
    print(f'\t[++] Preprocessed dataframe file created: {preprocessed_output_file_df}')
    df_preprocessed.to_csv(results_folder + os.path.sep + preprocessed_output_file_df, sep = ';', index = False)
      
    print('\t[++] Process file by HYGIEIA')
    
          
    call(["python", path_to_higieia + os.path.sep + exec_hygieia, results_folder, os.path.sep  + info["output_name"]])
    
    # hygieia config:         hygieia(path,model, opt_dup='2',opt_extra='2',verbose=1,interactive=False,
    #         sanitize=True, duplicates=True, inorganics=True, mixtures=True,
    #         salts=True, SEP=';', smiles_lab='SMILES',san_dups=True)
    
    
    processed_output_file_df = f'{info["output_name"]}.csv'
    print(f'\t[+++] Dataframe curated file created: {processed_output_file_df}')
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        