# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 18:45:54 2024

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

convert_ascii_data = True
data_folder =  '.' + os.path.sep + 'input_data'

results_folder = '.' + os.path.sep + 'results'

if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("Output directory created: ../results/")
else:
    print("Output directory already exists: ../results/")
    
clean_files_folder = '.' + os.path.sep + 'results' + os.path.sep + 'clean_files'

if not os.path.exists(clean_files_folder):
    os.makedirs(clean_files_folder)
    print("Output directory created: ../results/clean_files")
else:
    print("Output directory already exists: ../results/clean_files")


path_to_hygieia = r'F:\ProtoQSAR\GENERATE_MODELS_3_0\FullAsFunctions'
exec_hygieia = 'hygieia_mod.py'

config_file = 'info_original_datasets.xlsx'

##############################################################################

config_df = pd.read_excel('..' + os.path.sep + config_file)

config_df = config_df[config_df['proceed'] ==  'yes']

config_df_sel_cols = config_df[['dataset',
       'curation_process', 'sheet', 'SMILES_col', 'ID_col', 'NAME_col',
       'CAS_col', 'NAME_noiso', 'NAME_iso', 'CAS_noiso', 'CAS_iso', 'encoding', 'experimental_column', 'units',
       'output_name', 'require_additionalprocessing', 'proceed']]


dictios_dfs = config_df_sel_cols.set_index('dataset').T.to_dict()






for dataset, info in dictios_dfs.items():

    print(f'\n[+] Analysing "{dataset}" dataset')
    
    
    pathfile = data_folder + os.path.sep + dataset
    
    file_convertedname = dataset.split('.')[0]
    pathfile_converted = data_folder + os.path.sep + f'{file_convertedname}-converted.sdf'
    
    use_smiles_from_unique_source = True
    
   
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
        
        smiles_source = 'SMILES_from_mol'
            
            
    elif info['curation_process'] == 'curation2' or info['curation_process'] == 'curation5':
        
        if info['require_additionalprocessing'] == 'require_additionalprocessing_3':
        
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], header = 1)
            df.dropna(subset = [info['experimental_column']], inplace = True)
            df.replace(r'^>', '', regex=True, inplace = True)
            df[info['experimental_column']] = df[info['experimental_column']].astype(float)
        
        elif info['require_additionalprocessing'] == 'require_additionalprocessing_4':
            
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], header = 1)
            
        elif info['require_additionalprocessing'] == 'require_additionalprocessing_5':
            
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], engine='xlrd')
        
        else:
            df = pd.read_excel(pathfile, sheet_name = info['sheet'])
        
        original_output_file_df = f'{info["output_name"]}-original.csv'
        print(f'\t[++] Dataframe file created: {original_output_file_df}')
        df.to_csv(results_folder + os.path.sep + original_output_file_df, sep = ';', index = False)
        
        
        if info['require_additionalprocessing'] == 'require_additionalprocessing_1':
            df.rename(columns={'SMILES': 'treated_SMILES'}, inplace = True)
            
            smiles_source = info['SMILES_col']
            
        if info['require_additionalprocessing'] == 'require_additionalprocessing_2':
            allowed_endpoints = ['pKa', 'Acidic pKa', 'pKa(1)']
            df = df[df['Endpoint'].isin(allowed_endpoints)]
            df = df[df['Qualifier'] == 'No qualifier']
            
        if info['require_additionalprocessing'] == 'require_additionalprocessing_4':
            df = df[df['SDg'] < 2]
            
        smiles_source = info['SMILES_col']
            

    elif info['curation_process'] == 'curation3':
        path_test = data_folder + os.path.sep + dataset.split('|')[0] + '.sdf'
        df_test = PandasTools.LoadSDF(path_test)        
        path_train = data_folder + os.path.sep + dataset.split('|')[1] + '.sdf'
        df_train = PandasTools.LoadSDF(path_train)
        
        df = pd.concat([df_train,df_test], axis = 0)

        original_output_file_df = f'{info["output_name"]}-original.csv'
        print(f'\t[++] Dataframe file created: {original_output_file_df}')
        df.to_csv(results_folder + os.path.sep + original_output_file_df, sep = ';', index = False)  

                  
        if 'SMILES' in df.columns:
            df.rename(columns={'SMILES': 'orig_SMILES'}, inplace = True)
          
        
        smiles = [Chem.MolToSmiles(mol) for mol in df['ROMol']]
        df.insert(0,'SMILES',smiles)   
        
        smiles_source = 'SMILES_from_mol'

    elif info['curation_process'] == 'curation4':
        
        df = pd.read_csv(pathfile, sep = ',')
        
        smiles_source = info['SMILES_col']
        
    elif info['curation_process'] == 'curation6':
        
        df = pd.read_excel(pathfile, sheet_name = info['sheet'])
        
        use_smiles_from_unique_source = False
        
        sources = []
        
        selected_smiles = []
        
        smiles_source = []

        if pd.notnull(info['CAS_iso']):
            sources.append('IsofromCas')            
        if pd.notnull(info['NAME_iso']):
            sources.append('IsofromSmiles')        
        if pd.notnull(info['CAS_noiso']):
            sources.append('noIsofromCas') 
        if pd.notnull(info['NAME_noiso']):
            sources.append('noIsofromSmiles')    
        
        print(sources)
        
        for idx, row in df.iterrows():
            print(row)
            # if pd.notnull(row[info['CAS_iso']]):
            #     print('mimimimi')
            #     selected_smiles.append(row[info['CAS_iso']])
            #     smiles_source.append('isoSMILES_from_CAS')
            
            if pd.notnull(row[info['NAME_iso']]):
                selected_smiles.append(row[info['NAME_iso']])
                smiles_source.append('isoSMILES_from_NAME')
            
            # elif pd.notnull(row[info['CAS_noiso']]):
            #     selected_smiles.append(row[info['CAS_noiso']])
            #     smiles_source.append('noisoSMILES_from_CAS')
                
            elif pd.notnull(row[info['NAME_noiso']]):
                selected_smiles.append(row[info['NAME_noiso']])
                smiles_source.append('noisoSMILES_from_NAME')
            else:
                selected_smiles.append(np.nan)
                smiles_source.append(np.nan)
                
        # if info['require_additionalprocessing'] == 'require_additionalprocessing_6':
        #     df[info['experimental_column']] = df[info['experimental_column']]/100
            

        
        
#%%

    print(use_smiles_from_unique_source)

    if use_smiles_from_unique_source == True:    
        
        df.rename(columns = {info['SMILES_col']: 'SMILES'}, inplace = True)
        df.rename(columns = {info['NAME_col']: 'NAME'}, inplace = True)
    
        df.rename(columns = {info['CAS_col']: 'CAS'}, inplace = True)
        df.rename(columns = {info['ID_col']: 'orig_ID'}, inplace = True)
        
        df['smiles_source'] = [smiles_source]*df.shape[0]
        
    else:
        df.rename(columns = {'SMILES': 'SMILES_old'}, inplace = True)
        df['SMILES'] = selected_smiles
        df['smiles_source'] = smiles_source
        
        print('mami')
        # continue

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
        
    cols_to_maintain = ['ID', 'orig_ID', 'SMILES', info["experimental_column"], 'CAS', 'NAME', 'UNITS', 'smiles_source'] 
    
    df_preprocessed = df[cols_to_maintain]
    
    df_preprocessed.rename(columns={info["experimental_column"]: 'y'}, inplace = True)
   
    preprocessed_output_file_df = f'{info["output_name"]}-preprocessed.csv'
    print(f'\t[++] Preprocessed dataframe file created: {preprocessed_output_file_df}')
    df_preprocessed.to_csv(results_folder + os.path.sep + preprocessed_output_file_df, sep = ';', index = False)
      
    print('\t[++] Process file by HYGIEIA')
    
          
    call(["python", path_to_hygieia + os.path.sep + exec_hygieia, results_folder, os.path.sep  + info["output_name"]])
    
    # hygieia config:         hygieia(path,model, opt_dup='2',opt_extra='2',verbose=1,interactive=False,
    #         sanitize=True, duplicates=True, inorganics=True, mixtures=True,
    #         salts=True, SEP=';', smiles_lab='SMILES',san_dups=True)
    
    
    processed_output_file_df = f'{info["output_name"]}.csv'
    print(f'\t[+++] Dataframe curated file created: {processed_output_file_df}')
    
    df_final = pd.read_csv(results_folder + os.path.sep + processed_output_file_df,sep=';' )
    
    df_preprocessed.to_csv(clean_files_folder + os.path.sep + processed_output_file_df, sep = ';', index = False)
        
        
        
        
        
        
        
    
        
        
        
        
        
        
        
        