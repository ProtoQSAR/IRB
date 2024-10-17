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


# path_to_hygieia = '/home/carmenortiz/Escritorio/generate_models/FullAsFunctions/' #Carmen
path_to_hygieia = r'F:\ProtoQSAR\GENERATE_MODELS_3_0\FullAsFunctions' #Eva
exec_hygieia = 'hygieia_mod.py'

config_file = 'info_original_datasets_IRB.xlsx'
##############################################################################



config_df = pd.read_excel('..' + os.path.sep + config_file)

config_df = config_df[config_df['proceed'] ==  'yes']

config_df_sel_cols = config_df[['dataset',
       'curation_process', 'sheet', 'SMILES_col', 'ID_col', 'NAME_col',
       'CAS_col', 'NAME_noiso', 'NAME_iso', 'CAS_noiso', 'CAS_iso', 'encoding', 'experimental_column', 'units',
       'output_name', 'require_additionalprocessing', 'other_col', 'proceed']]
print(config_df)

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
            
        elif info['require_additionalprocessing'] == 'require_additionalprocessing_9':
            
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], skiprows=[1019, 1020, 1021], skipfooter=25)            

        elif info['require_additionalprocessing'] == 'require_additionalprocessing_10':
            
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], header = 5) 
            df['n_np'].replace({'p': 1, 'n': 0}, inplace = True)
             
        
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
        
        if info['require_additionalprocessing'] == 'require_additionalprocessing_8':
            df = pd.read_excel(pathfile, sheet_name = info['sheet'], header = 3)
        else:
            df = pd.read_excel(pathfile, sheet_name = info['sheet'])
        
        use_smiles_from_unique_source = False
        
        sources = []
        
        selected_smiles = []
        
        smiles_source = []

        
        
        
        for idx, row in df.iterrows():
            
            possible_smiles = []
            possible_sources = []
            
            if info['CAS_iso'] in row and pd.notnull(row[info['CAS_iso']]):
                
                
                try:
                    possible_smiles_san = Chem.MolToSmiles(Chem.MolFromSmiles(row[info['CAS_iso']]))
                except:
                    possible_smiles_san = np.nan

                possible_smiles.append(possible_smiles_san)
                possible_sources.append('isoSMILES_from_CAS')
            
            if info['NAME_iso'] in row and pd.notnull(row[info['NAME_iso']]):
                
                
                try:
                    possible_smiles_san = Chem.MolToSmiles(Chem.MolFromSmiles(row[info['NAME_iso']]))
                except:
                    possible_smiles_san = np.nan
                    
                possible_smiles.append(possible_smiles_san)
                possible_sources.append('isoSMILES_from_NAME')
            
            if info['CAS_noiso'] in row and pd.notnull(row[info['CAS_noiso']]):
                
                
                try:
                    possible_smiles_san = Chem.MolToSmiles(Chem.MolFromSmiles(row[info['CAS_noiso']]))
                except:
                    possible_smiles_san = np.nan
                    
                possible_smiles.append(possible_smiles_san)
                possible_sources.append('noisoSMILES_from_CAS')
                
            if info['NAME_noiso'] in row and pd.notnull(row[info['NAME_noiso']]):
                
                
                try:
                    possible_smiles_san = Chem.MolToSmiles(Chem.MolFromSmiles(row[info['NAME_noiso']]))
                except:
                    possible_smiles_san = np.nan
                possible_smiles.append(possible_smiles_san)
                possible_sources.append('noisoSMILES_from_NAME')
            
                
            
            
           
            
            selection = next(([smi, source] for smi, source in zip(possible_smiles, possible_sources) if pd.notnull(smi)), [np.nan, np.nan])
            
            
            selected_smiles.append(selection[0])
            smiles_source.append(selection[1])
      
            
                
        
                
        if info['require_additionalprocessing'] == 'require_additionalprocessing_6':
            df[info['experimental_column']] = df[info['experimental_column']]/10
            
        if info['require_additionalprocessing'] == 'require_additionalprocessing_7':
            
            df[info['experimental_column']] = pd.to_numeric(df[info['experimental_column']], errors='coerce')

        
    elif info['curation_process'] == 'curation7':
            
        df_orig = pd.read_excel(pathfile, sheet_name = info['sheet'], header = 1)
        
        df_orig.columns = ['id', 'compounds', 'melting point (°C)_ortho', 'melting point (°C)_meta', 'melting point (°C)_para',
       'SMILES_ortho', 'SMILES_meta', 'SMILES_para']
        
        df_orig[['melting point (°C)_ortho','melting point (°C)_meta', 'melting point (°C)_para']] = df_orig[['melting point (°C)_ortho','melting point (°C)_meta', 'melting point (°C)_para']].replace('c', np.nan)
        df_orig[['melting point (°C)_ortho','melting point (°C)_meta', 'melting point (°C)_para']] = df_orig[['melting point (°C)_ortho','melting point (°C)_meta', 'melting point (°C)_para']].replace('b', '', regex = True)
        df_orig[['melting point (°C)_ortho','melting point (°C)_meta', 'melting point (°C)_para']] = df_orig[['melting point (°C)_ortho','melting point (°C)_meta', 'melting point (°C)_para']].astype(float)
        
        
        
        df_orig.rename(columns = {info['NAME_col']: 'orig_NAME'}, inplace = True)
        
        
        into_dict = {}
        
        for idx, row in df_orig.iterrows():
            
            labels = ['o-', 'm-', 'p-']
            
            for label in labels:
                
                new_name = label + row['orig_NAME']
                into_dict[new_name] = {}
                
                if label == 'o-':
                    into_dict[new_name]['orig_ID'] = row[info['ID_col']]
                    into_dict[new_name]['SMILES'] = row['SMILES_ortho']
                    into_dict[new_name]['y'] = row['melting point (°C)_ortho']
                   
                elif label == 'm-':
                    into_dict[new_name]['orig_ID'] = row[info['ID_col']]
                    into_dict[new_name]['SMILES'] = row['SMILES_meta']
                    into_dict[new_name]['y'] = row['melting point (°C)_meta']
                    
                elif label == 'p-':
                    into_dict[new_name]['orig_ID'] = row[info['ID_col']]
                    into_dict[new_name]['SMILES'] = row['SMILES_para']
                    into_dict[new_name]['y'] = row['melting point (°C)_para']
                    

                    
        
        df = pd.DataFrame.from_dict(into_dict).T
        smiles_source = 'orig_SMILES_col'
        df.index.name = 'NAME'
        df.reset_index(inplace = True)
            
    elif info['curation_process'] == 'curation8':
        
       
        df = pd.read_csv(pathfile, sep = ';')
        
        if info['require_additionalprocessing'] == 'ignore_ID':
            df.rename(columns={'ID': 'old_ID'}, inplace = True)
        
        
        smiles_source = info['SMILES_col']
        
    if use_smiles_from_unique_source == True:    
        
        df.rename(columns = {info['SMILES_col']: 'SMILES'}, inplace = True)

        
        df['smiles_source'] = [smiles_source]*df.shape[0]
        
    else:
        df.rename(columns = {'SMILES': 'SMILES_old'}, inplace = True)
        df['SMILES'] = selected_smiles
        df['smiles_source'] = smiles_source
        
    
    if pd.notnull(info['other_col']):
        print(dataset)
        new_col = info['other_col'].split('_')[0]
        col_to_get = info['other_col'].split('_')[1]
        df[new_col] = df[col_to_get]
                       
    
    

    df.rename(columns = {info['NAME_col']: 'NAME'}, inplace = True)

    df.rename(columns = {info['CAS_col']: 'CAS'}, inplace = True)
    df.rename(columns = {info['ID_col']: 'orig_ID'}, inplace = True)



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
        df['NAME'] = df['NAME'].astype(str).str.replace('|', '_', regex = True)
        df['NAME'] = df['NAME'].astype(str).str.replace(';', ',', regex = True)
        
        
    if pd.notnull(info['other_col']):
                
        cols_to_maintain = ['ID', 'orig_ID', 'SMILES', info["experimental_column"], 'CAS', 'NAME', 'UNITS', 'smiles_source', new_col]
        
    else:
        
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
    
    
    processed_output_file_df = f'{info["output_name"]}'
    print(f'\t[+++] Dataframe curated file created: {processed_output_file_df}')
    
    df_final = pd.read_csv(results_folder + os.path.sep + processed_output_file_df + '.csv',sep=';' )
    
    df_final = df_final[df_final['y'].notna()]
    
    if pd.notnull(info['other_col']):
         cols_to_maintain2 = ['ID', 'orig_ID', 'SMILES', 'y', 'CAS', 'NAME', 'UNITS', 'smiles_source', new_col] 
    else:
        
        cols_to_maintain2 = ['ID', 'orig_ID', 'SMILES', 'y', 'CAS', 'NAME', 'UNITS', 'smiles_source'] 
    
    
    df_final = df_final[cols_to_maintain2]
    
    df_final.replace({'[nan, nan]':np.nan, '[nan, nan, nan]':np.nan, '[nan, nan, nan, nan]':np.nan}, inplace = True)
    
    df_final.to_csv(clean_files_folder + os.path.sep + processed_output_file_df + '_final.csv', sep = ';', index = False)
        
        
        
        
        
        
        
    
        
        
        
        
        
        
        
        
