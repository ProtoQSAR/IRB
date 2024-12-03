# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 10:14:07 2024

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
##############################################################################
################################ INITIAL VARIABLES ###########################

parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')
data_folder =  '..' + os.path.sep + 'compare_ONGOING' + os.path.sep + 'results'
config_file_files = 'config_predictions_files.xlsx'
config_file_tools = 'config_predictions_tools.xlsx'

##############################################################################

###############################  FUNCTIONS ###################################

def read_file(file, separator):

    pathfile = data_folder + os.path.sep + file

    extension = file[-3:]

    if extension == 'csv':
        df = pd.read_csv(pathfile, sep = separator)


    return df

def create_parts(df_output, config_df_tools, tool_to_manage, dataset):

    if tool_to_manage == 'SwissADME_tool':

        df_output = df_output[['SMILES','NAME']]
        to_drop = []
        for idx, row in df_output.iterrows():
            if len(row['SMILES']) >= 196: #(considering numbering)
                to_drop.append(idx)

        if len(to_drop) !=0:

            df_output.drop(to_drop, inplace = True)
            print(f'\t\tsome molecules [{len(to_drop)}] have been dropped (only 200 characters allowed)')

    # if tool_to_manage == 'ADMETLab3_tool':

    #     to_drop = []
    #     for idx, row in df_output.iterrows():
    #         # print(row['SMILES'], len(row['SMILES']))
    #         if len(row['SMILES']) >= 200:
    #             to_drop.append(idx)

    #     if len(to_drop) >1:
    #         print(f'\t\tsome molecules [{len(to_drop)}] have been dropped (only 200 characters allowed)')
    #         df_output.drop(to_drop, inplace = True)



    #split dataframe
    max_rows = int(config_df_tools.loc[tool_to_manage,'Threshold    (# molecules)'])
    part_dataframes = []
    while len(df_output) > max_rows:
        top = df_output[:max_rows]
        part_dataframes.append(top)
        df_output = df_output[max_rows:]
    else:
        part_dataframes.append(df_output)

    #create output dirs
    inputs_to_predict_folder = results_folder + os.path.sep + tool_to_manage + os.path.sep + 'inputs'

    create_dir_for_io(inputs_to_predict_folder)


    #save parts

    separator_specified = config_df_tools.loc[tool_to_manage,'separator']

    if pd.notnull(separator_specified):

        if separator_specified == 'space':

            separator_for_csv = ' '
        else:
            separator_for_csv = separator_specified



    else:
        separator_for_csv = ','

    print('\t\t[+++] The following files have been created:')

    for number, part_dataframe in enumerate(part_dataframes):

        tag_dataset = ''.join(dataset.split('_'))
        tag_tool = tool_to_manage.split('_')[0]

        file_name = f'{tag_dataset}_input{tag_tool}_part{number}.csv'
        part_dataframe.to_csv(inputs_to_predict_folder + os.path.sep + file_name, sep=separator_for_csv, index=False)
        # print(part_dataframe.shape)
        print(f'\t\t\t{file_name} __ num mols: {part_dataframe.shape[0]}')

    outputs_to_predict_folder = results_folder + os.path.sep + tool_to_manage + os.path.sep + 'predictions'

    create_dir_for_io(outputs_to_predict_folder)



def create_dir_for_io(folder):

    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"\t\tDirectory created: {folder}")
    else:
        print(f"\t\tDirectory already exists: {folder}")



def prepareqsartoolbox(row):
    row = 'MoleculeID' + str(row)
    
    return row

def elimieisomery(row):
    
    from rdkit import Chem
    
    noisosmiles = Chem.MolToSmiles(Chem.MolFromSmiles(row), isomericSmiles=False)
    
    return noisosmiles


##############################################################################

############################### CONFIGURATION FILES ##########################

results_folder = '.' + os.path.sep + 'prediction_files'

create_dir_for_io(results_folder)

# for tool parameters
config_df_tools = pd.read_excel(config_file_tools, sheet_name = 'tools_configuration', header = 1)

config_df_tools = config_df_tools[['Tool name', 'Tool_column', 'Threshold    (# molecules)', 'configured', 'separator']]

config_df_tools.set_index('Tool_column', inplace = True)

# for tools enpoints

df_tools_enpoints = pd.read_excel(config_file_tools, sheet_name = 'properties_to_prepare', header = 1)

df_tools_enpoints = df_tools_enpoints.drop('Category', axis = 1)

dict_tools_enpoints = df_tools_enpoints.set_index('Tool_column').T.to_dict(orient = 'list')


# for input files
config_df = pd.read_excel(config_file_files, sheet_name = 'input_files')

config_df = config_df[config_df['proceed_create_input_files'] ==  'yes']

config_df.set_index('dataset_ID', inplace = True)

tool_cols = [tool for tool in config_df.columns if 'tool' in tool]


#%%
##############################################################################

for dataset, info in config_df.iterrows():

    print(f'\n[+] Analysing "{dataset}" dataset')

    endpoint = f"{dataset.split('_')[0]}_{dataset.split('_')[1]}"


    file = info['file']
    separator = info['separator']

    df = read_file(file, separator)

    tools_to_manage = [tool for tool in list(config_df_tools.index)]

    for tool_to_manage in tools_to_manage:

        if config_df_tools.loc[tool_to_manage,'configured'] == 'yes':

            if tool_to_manage in dict_tools_enpoints[endpoint]:


                print(f"\t[++] Creating files for {tool_to_manage}")

                if tool_to_manage == 'ProtoPRED_tool': # compulsory columns SMILES and EC number (for ID)

                    try:
                        df_output = df[[info['Other Regulatory ID'], info['SMILES']]]
                    except:
                        print(f"\tWARNING!!!! 'SMILES' and 'EC number' (with ID) columns are compulsory for {tool_to_manage}")

                    df_output.rename(columns = {'ID_all' : 'Other Regulatory ID'}, inplace = True)

                    if pd.notnull(info['CAS']):
                        df_output.loc[:,'CAS'] = list(df[info['CAS']])
                    if pd.notnull(info['Structural formula']):
                        df_output.loc[:,'Structural formula'] = list(df[info['Structural formula']])
                    if pd.notnull(info['Chemical name']):
                        df_output.loc[:,'Chemical name'] = list(df[info['Chemical name']])


                    create_parts(df_output, config_df_tools, tool_to_manage, dataset)

                if tool_to_manage == 'ADMETLab3_tool' or tool_to_manage == 'PKCSM_tool': # compulsory columns SMILES

                    try:
                        df_output = df[[info['SMILES']]]
                    except:
                        print(f"\tWARNING!!!! 'SMILES' column is compulsory for {tool_to_manage}")


                    create_parts(df_output, config_df_tools, tool_to_manage, dataset)


                if tool_to_manage == 'SwissADME_tool':

                    try:
                        df_output = df[[info['NAME'], info['SMILES']]]
                    except:
                        print(f"\tWARNING!!!! 'NAME' and 'SMILES' columns are compulsory for {tool_to_manage}")

                    df_output.rename(columns = {'ID_all' : 'NAME'}, inplace = True)


                    create_parts(df_output, config_df_tools, tool_to_manage, dataset)



                if tool_to_manage == 'vNNADMET_tool': # compulsory columns SMILES and NAME (for ID)

                    try:
                        df_output = df[[info['NAME'], info['SMILES']]]
                    except:
                        print(f"\tWARNING!!!! 'NAME' and 'SMILES' columns are compulsory for {tool_to_manage}")


                    df_output.rename(columns = {'ID_all' : 'NAME'}, inplace = True)

                        #create output dirs
                    inputs_to_predict_folder = results_folder + os.path.sep + tool_to_manage + os.path.sep + 'inputs'

                    create_dir_for_io(inputs_to_predict_folder)

                    outputs_to_predict_folder = results_folder + os.path.sep + tool_to_manage + os.path.sep + 'predictions'

                    create_dir_for_io(outputs_to_predict_folder)

                    tag_dataset = ''.join(dataset.split('_'))
                    tag_tool = tool_to_manage.split('_')[0]

                    file_name = f'{tag_dataset}_input{tag_tool}.csv'

                    df_output.to_csv(inputs_to_predict_folder + os.path.sep + file_name, index=False)

                    print(f'\t\t\t{file_name}')
                    
                if tool_to_manage == 'QSARToolbox_tool':
                    

                    df_output = df[[info['NAME'], info['SMILES']]]
                    
                

                    
                    
                    df_output['NAME'] = df_output.apply(lambda row: prepareqsartoolbox(row[info['NAME']]), axis=1)
                    
                    df_output['nonisoSMILES'] = df_output.apply(lambda row: elimieisomery(row[info['SMILES']]), axis=1)
                    
                    df_output = df_output[['NAME', 'nonisoSMILES']]
                    
                    df_output.columns = ['NAME', 'SMILES']
                    
                    inputs_to_predict_folder = results_folder + os.path.sep + tool_to_manage + os.path.sep + 'inputs'

                    create_dir_for_io(inputs_to_predict_folder)

                    tag_dataset = ''.join(dataset.split('_'))
                    tag_tool = tool_to_manage.split('_')[0]
                    
                    file_name = f'{tag_dataset}_input{tag_tool}.txt'
                    
                    df_output.to_csv(inputs_to_predict_folder + os.path.sep + file_name, index=False, sep ='\t')

                    print(f'\t\t\t{file_name}') 


                    rows_red = [14904, 15129,12137,12280,6899,6902,4817,4829,6641,6642,6608,6613]
                    
                    names = ['MoleculeID'+str(x) for x in rows_red]
                    
                    df_reduced = df_output[df_output['NAME'].isin(names)]


                    
                    file_name = f'{tag_dataset}_input{tag_tool}_reduced.txt'
                    
                    
                    
                    
                    df_reduced.to_csv(inputs_to_predict_folder + os.path.sep + file_name, index=False, sep ='\t')
                    

                    

                    
                    

        #     else:
        #         print(f"\t[++] {tool_to_manage} do not predict the endpoint {endpoint}. No files have been created")




        # else:
        #     print(f"\t[++] {tool_to_manage} need to be configured")
