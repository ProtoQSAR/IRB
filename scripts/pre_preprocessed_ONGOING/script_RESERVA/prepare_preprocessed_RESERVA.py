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
    


##############################################################################

dictios_dfs = {
    # 'EPI_AOP_Data_SDF.sdf' : {},
    # 'EPI_BCF_Data_SDF.sdf': {},
    # 'EPI_BioHC_Data_SDF.sdf': {},
    # 'EPI_Biowin_Data_SDF.sdf': {},
    'EPI_Boil_Pt_Data_SDF.sdf': {'encoding': 'ISO-8859-1', 'experimental_column' : 'ExpBP','units' : 'C' , 'output_name' : 'PhyspropBP'},
    # 'EPI_Henry_Data_SDF.sdf': {},
    # 'EPI_KM_Data_SDF.sdf': {},
    # 'EPI_KOA_Data_SDF.sdf': {},
    'EPI_Kowwin_Data_SDF.sdf': {'encoding': 'ISO-8859-1', 'experimental_column' : 'Kow','units' : 'logKow' , 'output_name' : 'PhyspropLogP'},
    'EPI_Melt_Pt_Data_SDF.sdf': {'encoding': 'ISO-8859-1', 'experimental_column' : 'MP','units' : 'C' , 'output_name' : 'PhyspropMP'},
    # 'EPI_PCKOC_Data_SDF.sdf': {},
    'EPI_VP_Data_SDF.sdf': {'encoding': 'utf8', 'experimental_column' : 'VP','units' : 'mmHg' , 'output_name' : 'PhyspropVP'},
    # 'EPI_WaterFrag_Data_SDF.sdf': {},
    'EPI_Wskowwin_Data_SDF.sdf': {'encoding': 'ISO-8859-1', 'experimental_column' : 'WS','units' : 'mg/L' , 'output_name' : 'PhyspropWS'},
    }

for dataset, info in dictios_dfs.items():

    print(f'\n[+] Analysing "{dataset}" dataset')
    
    extension = dataset.split('.')[1]
    
    pathfile = data_folder + os.path.sep + dataset
    
    file_convertedname = dataset.split('.')[0]
    pathfile_converted = data_folder + os.path.sep + f'{file_convertedname}_converted.sdf'
    
   
    if extension == 'sdf':   
        # continue
       
       
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

        if 'ID' in df.columns:
            df.rename(columns={'ID': 'orig_ID'}, inplace = True)
        
        smiles = [Chem.MolToSmiles(mol) for mol in df['ROMol']]
        df.insert(0,'SMILES',smiles)
        
        df.index.name = 'ID'
        df.reset_index(inplace = True)
        
        df.insert(df.shape[1],'experimental_column',info["experimental_column"])
        df.insert(df.shape[1],'units',info["units"])
        
        entire_output_file_df = f'{info["output_name"]}-entiredata.csv'
        print(f'\t[++] Dataframe file created: {entire_output_file_df}')
        df.to_csv(results_folder + os.path.sep + entire_output_file_df, sep = ';', index = False)
        
        cols_to_maintain = ['ID', 'SMILES', info["experimental_column"]]
        
        if 'CAS' in list(df.columns):
            cols_to_maintain.append('CAS')
        if 'NAME' in list(df.columns):
            cols_to_maintain.append('NAME')
        
        df_preprocessed = df[cols_to_maintain]
        
        df_preprocessed['NAME'] = df_preprocessed['NAME'].str.replace('|', '_', regex = True)
        df_preprocessed['NAME'] = df_preprocessed['NAME'].str.replace(';', ',', regex = True)
        
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        