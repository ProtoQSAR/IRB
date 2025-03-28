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
import re   #"""importar regular expresions"""
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
config_file_files = 'config_predictions_files.xlsx'
config_file_tools = 'config_predictions_tools.xlsx'

delimite_endpoints = True
merge_SITOLUB=True
merge_ONTOX = False #delimite_endpoints should be True

##############################################################################

###############################  FUNCTIONS ###################################


def eliminar_version(model):

    return re.sub(r"\(version \d*\.\d*\.\d*\)?", "", model).strip()
##"""Modified by Vicente 09/12/2024. skiprow (detects header of VEGA atomatically)"""
def detectar_skiprow(file_path, header_start="No.",return_models=False):
    models = []
    with open(file_path, "r", encoding="ANSI") as file:
        lines = file.readlines()


    for i, line in enumerate(lines):
        if header_start in line:
            #  print(i)
            if return_models:
                models = [lines[j].strip() for j in range(1, i-2)]
                models = [eliminar_version(model) for model in models]

                #  print(models)
                return models
            return i
    return None if not return_models else models

def create_dir_for_io(folder):

    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"\t\tDirectory created: {folder}")
    else:
        print(f"\t\tDirectory already exists: {folder}")


def collect_files(folder):
    for root, dirs, files_folder in os.walk(folder):
            continue

    files_by_property = [x.split('.')[0] for x in files_folder if x.startswith(complete_endpoint) ]

    return files_by_property

# def read_file_predictions(file, path_for_files, equivalence_in_tool):

#     if tool == 'ADMETLab3_tool':
#         df = pd.read_csv(path_for_files + os.path.sep + file + '.csv')

#     if tool == 'ProtoPRED_tool':

#     return df


def merge_files(list_files, path_for_files, tool, equivalence_in_tool, mode):
    
   

    list_to_merge = []
    
    models = []
    for file in list_files:

        if tool == 'ADMETLab3_tool':
            df = pd.read_csv(path_for_files + os.path.sep + file + '.csv')

            if 'molstr' in df.columns: # this generates conflicts and is the mol representation
                df = df.drop('molstr', axis = 1)
            list_to_merge.append(df)
            

        elif tool == 'ProtoPRED_tool':

            if mode == 'inputfiles':
                df = pd.read_csv(path_for_files + os.path.sep + file + '.csv', sep = ';')
            elif mode == 'outputfiles':

                df = pd.read_excel(path_for_files + os.path.sep + file + '.xlsx', sheet_name= equivalence_in_tool)

            list_to_merge.append(df)

        elif 'SwissADME' in tool:

            df = pd.read_csv(path_for_files + os.path.sep + file + '.csv', sep = ',')

            list_to_merge.append(df)

        elif 'VEGA' in tool: #"""Modified by Vicente 09/12/2024. Add Vega tool"""
            
            skiprows = detectar_skiprow(path_for_files + os.path.sep + file + '.txt')
            models = detectar_skiprow(path_for_files + os.path.sep + file + '.txt', return_models=True)


            if mode=="inputfiles":
                df = pd.read_csv(
                    path_for_files + os.path.sep + file + '.txt',
                    sep="\t",
                    encoding="ANSI",
                    skiprows=skiprows,
                    header=None
                )
                df["SMILES"]=df
                df=df[["SMILES"]]

            else:
                df = pd.read_csv(
                    path_for_files + os.path.sep + file + '.txt',
                    sep="\t",
                    encoding="ANSI",
                    skiprows=skiprows,
                )
                rename_columns = {
                    df.columns[i]: f'Pred_{df.columns[i]}' for i in range(4, len(df.columns), 2)
                }
                df.rename(columns=rename_columns, inplace=True)
                new_columns_exp = {}
                new_columns_AD = {}


                for i in range(3, len(df.columns), 2):
                    col_name = df.columns[i]


                    new_columns_exp[f'Exp_{col_name}'] = np.where(
                    df[col_name].str.contains("EXPERIMENTAL"), 
                    df[col_name].str.extract(r'(-?[\d.]+)')[0].astype(float).apply(lambda x: round(x, 2)), 
                    np.nan  
                    )
                    # new_columns_AD[f'AD_{col_name}'] =  df[col_name].str.extract(r'\((.*?)\)') 
                    new_columns_AD[f'AD_{col_name}'] = np.where(
                        df[col_name].str.contains(r'\(LOW reliability\)'),  
                        0,  
                        1   
                    )

                # Añadir las nuevas columnas al DataFrame
                for col_name, col_data in new_columns_exp.items():
                    df[col_name] = col_data

                for col_name, col_data in new_columns_AD.items():
                    df[col_name] = col_data
                filtered_columns = df.filter(regex=r'^(tId|SMILES|Exp|AD|Pred)', axis=1)
                df = filtered_columns
                # print(df.head())
                # print(df.columns)

            list_to_merge.append(df)
            
        elif 'ECOSAR' in tool:
            df=pd.read_csv(path_for_files + os.path.sep + file + '.csv', sep = ";", )
            
            list_to_merge.append(df)
            
        elif tool == 'EXPERIMENTAL_tool': # 'EXPERIMENTAL' in tool: #"""Modified by Vicente 12/12/2024. Add Experimental tool"""
            """_summary_ añadir valores experimentales, de clasificación o regresión + columna SMILES """
            df = pd.read_csv(path_for_files + os.path.sep + file + '.csv', sep = ";", )

            list_to_merge.append(df)
            
        elif tool == 'PKCSM_tool':

            df = pd.read_csv(path_for_files + os.path.sep + file + '.csv', sep = '\t', lineterminator='\n')

            if 'SMILES' in list(df.columns):

                df['SMILES'].replace('\r','', inplace= True, regex = True)


            list_to_merge.append(df)
            
    

    merged_file = pd.concat(list_to_merge)
    # print(merged_file,"PATATA")
    merged_file.reset_index(inplace=True)

    file_name = list_files[0].split('_part')[0]
    
    return merged_file, file_name, models



def conversion(row, units, desired_units):


    if units == 'proba_BBBnormal' and desired_units == 'class_BBB':
        if row > 0.5: return 1
        else: return 0

    elif units == 'proba_Psubnormal' and desired_units == 'class_Psub':
        if row > 0.5: return 1
        else: return 0

    elif units == 'proba_Pinhnormal' and desired_units == 'class_Pinh':
        if row > 0.5: return 1
        else: return 0


    elif units == 'logBBB' and desired_units == 'class_BBB':

        if row >= (-1): return 1
        else: return 0

    elif units == 'mmHg' and desired_units == 'LogmmHg':

        return np.log10(row)

    elif units == 'proba_F30inverse' and desired_units == 'class_F30':

        if row > 0.5: return 0
        else: return 1


    elif units == 'ratioF10' and desired_units == 'class_F30':

        if row >= 0.5: return 1
        else: return 0


    elif units == 'proba_HIAinverse' and desired_units == 'class_HIA':

        if row > 0.5: return 0
        else: return 1


    elif units == 'lowhighclass_HIA' and desired_units == 'class_HIA':

        if row == 'High': return 1
        else: return 0





    elif units == '%FU' and desired_units == 'ratioFU':

        return row/100




    elif units == 'M' and desired_units == 'logM':

        return np.log(row)

    elif units == 'logKp(cm/s)' and desired_units ==    'logKp(cm/h)':

        antilog = 10**row

        converted = antilog*3600

        convertedlog = np.log10(converted)

        return convertedlog

    elif units == 'log10-6(cm/s)' and desired_units == 'log(cm/s)':

        value = 10**row

        value2 = value*10**(-6)

        res = np.log10(value2)

        return res





    else:

        if isinstance(row, str):

            row = row.replace('Yes', '1')
            row = row.replace('No', '0')


        return row

def ad_vnnadmet(row,equivalence_in_tool):

    ads = 0 if row[f'{equivalence_in_tool}_restricted'] == 'NoPrediction' else 1


    return ads

def ad_ProtoPRED(row):

    ads = 0 if 'Outside' in row else 1

    return ads

def ad_ADMETlab(row, equivalence_in_tool, desired_units):



    proba = row[equivalence_in_tool]


    if proba >= 0.7 or proba <= 0.3:
        ad_ADMETlab3 = 1
    else:
        ad_ADMETlab3 = 0

    return ad_ADMETlab3



def ts_ProtoPRED(row):

    ts = 0 if (row == '-' or pd.isnull(row) or row == 'N/A') else 1
    return ts
def exp_ProtoPRED():
    np.where(
                    df["Experimental numerical"].str.contains("-"), 
                    np.nan, 
                    df["Experimental numerical"].astype(float))  
                    

def sanitize_smiles(row):
    mol = Chem.MolFromSmiles(row)


    san_smi = Chem.MolToSmiles(mol)

    return san_smi


##############################################################################

############################### CONFIGURATION FILES ##########################

results_folder = '.' + os.path.sep + 'prediction_files'

# for tool parameters
config_df_tools = pd.read_excel(config_file_tools, sheet_name = 'tools_configuration', header = 1)

config_df_tools = config_df_tools[['Tool_column', 'smiles_output', 'id_output']]

config_df_tools.set_index('Tool_column', inplace = True)


# for tools endpoints

df_tools_enpoints = pd.read_excel(config_file_tools, sheet_name = 'properties_to_proccess', header = 1)

df_tools_enpoints = df_tools_enpoints.drop('Category', axis = 1)

dict_tools_enpoints = df_tools_enpoints.set_index('Tool_column').T.to_dict(orient = 'list')



# for tools equivalence in names

df_tools_enpoint_names = pd.read_excel(config_file_tools, sheet_name = 'properties_by_tools_equivalence', header = 1)

df_tools_enpoint_names = df_tools_enpoint_names.drop('Category', axis = 1)

df_tools_enpoint_names.set_index('Tool_column', inplace = True)

# for tools units

df_tools_enpoint_units = pd.read_excel(config_file_tools, sheet_name = 'properties_by_tools_units', header = 1)

df_tools_enpoint_units = df_tools_enpoint_units.drop('Category', axis = 1)

df_tools_enpoint_units.set_index('Tool_column', inplace = True)




# for input files
config_df = pd.read_excel(config_file_files, sheet_name = 'input_files')

config_df = config_df[config_df['proccess_predictions'] ==  'yes']

config_df.set_index('dataset_ID', inplace = True)



if merge_ONTOX == True:

    # for original_files to merge
    initial_data_folder =  '..' + os.path.sep + 'compare_ONGOING' + os.path.sep + 'results'

if merge_SITOLUB == True:

    # for original_files to merge
    initial_data_folder =  '..' + os.path.sep + 'compare_ONGOING' + os.path.sep + 'results'



##############################################################################






for dataset, info in config_df.iterrows():

    print(f'\n[+] Analysing "{dataset}" dataset')

    split_dataset_name = dataset.split('_')
    endpoint = split_dataset_name[1]
    complete_endpoint = ''.join(split_dataset_name[0:2])
    tag = '_'.join(split_dataset_name[0:2])

    # fist check step
    tools_to_process = [tool for tool in dict_tools_enpoints[tag] if pd.notnull(tool)]

    
    processed_dfs = {}
    validate_df={}
    val = pd.DataFrame()
    print('  [++] Merging prediction files')

    for tool in tools_to_process:
        if tool != 'VEGA_tool' and tool != "EXPERIMENTAL_tool" and tool != "ECOSAR_tool":

            print('melon')
            print(f'\t[+++] Analysing "{tool}" tool')

            tag_tool = tool.split('_')[0]

            if 'SwissADME' in tool:


                tool_folder = results_folder + os.path.sep + 'SwissADME_tool'

            else:

                tool_folder = results_folder + os.path.sep + tool


            files_to_predict_folder = tool_folder + os.path.sep + 'inputs'
            files_predicted_folder = tool_folder + os.path.sep + 'predictions'

            equivalence_in_tool = df_tools_enpoint_names.loc[tag,tool]
            # print(equivalence_in_tool,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")




            #collect files to be predicted
            files_topredict_property = collect_files(files_to_predict_folder)



            #collect files predicted
            files_predicted_property = collect_files(files_predicted_folder)

            if tool != 'vNNADMET_tool':
                ########### compare files ---> FIRST CHECK
                input_to_output = [name.replace('input','output') for name in files_topredict_property]
                different = set(input_to_output).difference(set(files_predicted_property))

                if len(different) > 0:

                    print(f'\t\t {len(files_topredict_property)} files in input')
                    print(f'\t\t {len(files_predicted_property)} files will be processed')

                    print('\t\t WARNING!!!! some input files have not prediction')
                    print('\t\t processing cannot continue')
                    break
                else:
                    print(f'\t\t {len(files_topredict_property)} files in input')
                    print(f'\t\t {len(files_predicted_property)} files will be processed')



                # merge input files

                merged_input, file_name_input, models = merge_files(files_topredict_property, files_to_predict_folder, tool, equivalence_in_tool, mode = 'inputfiles')
                # merge output files

                merged_output, file_name_output, models = merge_files(files_predicted_property, files_predicted_folder, tool, equivalence_in_tool, mode = 'outputfiles')

           

                smile_column_from_tool = config_df_tools.loc[tool, 'smiles_output']
                id_column_from_tool = config_df_tools.loc[tool, 'id_output']


                merged_output.rename(columns = {smile_column_from_tool : f'SMILES_{tag_tool}',id_column_from_tool : f'ID_{tag_tool}' }, inplace = True)

                # if pd.isnull(id_column_from_tool):
                #     merged_output.index.name = f'ID_{tag_tool}'
                #     merged_output.reset_index(inplace = True)

                ########### compare number of molecules ---> SECOND CHECK

                if merged_input.shape[0] != merged_output.shape[0]:
                    print(f'\t\t WARNING!!!! number of input molecules {merged_input.shape[0]} != number of output molecules {merged_output.shape[0]}')

                else:
                    print(f'\t\t {merged_output.shape[0]} molecules found')

            

            else:

                try:
                    restricted = [x for x in files_predicted_property if '_restricted' in x]
                    restricted_df = pd.read_csv(files_predicted_folder + os.path.sep + restricted[0] + '.csv')
                    restricted_df_columns = [f'{col}_restricted' for col in restricted_df]
                    restricted_df.columns = restricted_df_columns
                    restricted_df.set_index('Query_restricted')

                    unrestricted = [x for x in files_predicted_property if '_unrestricted' in x]
                    unrestricted_df = pd.read_csv(files_predicted_folder + os.path.sep + unrestricted[0] + '.csv')
                    unrestricted_df_columns = [f'{col}_unrestricted' for col in unrestricted_df]
                    unrestricted_df.columns = unrestricted_df_columns
                    unrestricted_df.set_index('Query_unrestricted')

                    merged_output = pd.concat([restricted_df, unrestricted_df], axis = 1)

                    smile_column_from_tool = config_df_tools.loc[tool, 'smiles_output']
                    id_column_from_tool = config_df_tools.loc[tool, 'id_output']

                    merged_output.rename(columns = {smile_column_from_tool : f'SMILES_{tag_tool}',id_column_from_tool : f'ID_{tag_tool}' }, inplace = True)


                    file_name_output = restricted[0].strip('_restricted')
                except:
                    print('\t\t WARNING!!!! some input files have not prediction')
                    print('\t\t processing cannot continue')
                    break


            merged_files_folder = results_folder + os.path.sep + tool + os.path.sep + 'merged_predictions_all'
            create_dir_for_io(merged_files_folder)
            merged_output.to_csv(merged_files_folder + os.path.sep + file_name_output + '_allcolumns.csv', sep = ';', index = False)



            if delimite_endpoints == True:






                if tool == 'ProtoPRED_tool':

                    ''' just an exception on WS to get the value in logM to avoid internal conversion'''
                    ###############################################################

                    if endpoint == 'WS':

                        # cols = ['Other Regulatory ID',	'SMILES',	'Experimental numerical (common units)',	'Predicted numerical (model units)','Applicability domain**'] """Modified by Vicente 09/12/2024. Motive: Change in column name of Protopred predictions"""
                        cols = ['Other Regulatory ID',	'SMILES',	'Experimental numerical (model units)',	'Predicted numerical (model units)','Applicability domain**']
                        merged_output.rename(columns = {'Other Regulatory ID' : f'ID_{tag_tool}',
                                                    'SMILES' : f'SMILES_{tag_tool}',
                                                    'Experimental numerical (common units)': f'exp_{tag_tool}_{endpoint}',
                                                    'Predicted numerical (model units)': f'pred_{tag_tool}_{endpoint}',

                                                    }, inplace = True)
                        # print(merged_output.columns)
                    else:

                        cols = ['Other Regulatory ID',	'SMILES',	'Experimental numerical',	'Predicted numerical','Applicability domain**']





                        merged_output.rename(columns = {'Other Regulatory ID' : f'ID_{tag_tool}',
                                                    'SMILES' : f'SMILES_{tag_tool}',
                                                    'Experimental numerical': f'exp_{tag_tool}_{endpoint}',
                                                    'Predicted numerical': f'pred_{tag_tool}_{endpoint}',

                                                    }, inplace = True)


                    ###############################################################
                    merged_output[f'AD_{tag_tool}_{endpoint}'] = merged_output.apply(lambda row: ad_ProtoPRED(row['Applicability domain**']), axis=1)
                    merged_output[f'in_TS_{tag_tool}_{endpoint}'] = merged_output.apply(lambda row: ts_ProtoPRED(row[f'exp_{tag_tool}_{endpoint}']), axis=1)


                    cols = [f'ID_{tag_tool}', f'SMILES_{tag_tool}', f'exp_{tag_tool}_{endpoint}',f'pred_{tag_tool}_{endpoint}',f'AD_{tag_tool}_{endpoint}',f'in_TS_{tag_tool}_{endpoint}' ]








                else:
                    if tool == 'vNNADMET_tool':


                        merged_output[f'AD_{tag_tool}_{endpoint}'] = merged_output.apply(lambda row: ad_vnnadmet(row, equivalence_in_tool), axis=1)

                        merged_output.rename(columns= {f'{equivalence_in_tool}_unrestricted' : f'pred_{tag_tool}_{endpoint}',
                                                    f'{equivalence_in_tool}_restricted' : f'pred_{tag_tool}_{endpoint}_forad'}, inplace =True)



                        cols = [f'ID_{tag_tool}', f'SMILES_{tag_tool}', f'pred_{tag_tool}_{endpoint}', f'pred_{tag_tool}_{endpoint}_forad', f'AD_{tag_tool}_{endpoint}' ]


                    else:

                        if pd.notnull(id_column_from_tool):

                            cols = [f'ID_{tag_tool}', f'SMILES_{tag_tool}', equivalence_in_tool]

                        

                        else:
                            cols = [f'SMILES_{tag_tool}', equivalence_in_tool]

                merged_output_selcols = merged_output[cols]


                merged_output_selcols.rename(columns = {equivalence_in_tool : f'pred_{tag_tool}_{endpoint}'}, inplace = True)



                predicted_units = df_tools_enpoint_units.loc[tag,tool]

                desired_units = config_df.loc[dataset,'desired_units']

                if tool == 'ADMETLab3_tool':
                    predicted_units = predicted_units.split(':')[0]

                    if 'class' in desired_units:

                        merged_output_selcols[f'AD_{tag_tool}_{endpoint}'] = merged_output.apply(lambda row: ad_ADMETlab(row, equivalence_in_tool, desired_units), axis=1)



                # print(desired_units, 'desired_units')

                # print(predicted_units, 'predicted_units')

                merged_output_selcols[f'units_{tag_tool}_{endpoint}'] = predicted_units

                merged_output_selcols[f'pred_{tag_tool}_{endpoint}_converted'] = merged_output_selcols.apply(lambda row: conversion(row[f'pred_{tag_tool}_{endpoint}'], predicted_units, desired_units), axis=1)

                merged_output_selcols[f'pred_{tag_tool}_{endpoint}_convertedunits'] = desired_units


                merged_files_folder_selcols = results_folder + os.path.sep + tool + os.path.sep + 'merged_predictions_selcols'
                create_dir_for_io(merged_files_folder_selcols)
                merged_output_selcols.to_csv(merged_files_folder_selcols + os.path.sep + file_name_output + f'_{complete_endpoint}.csv', sep = ';', index = False)

                merged_files_folder_selcols_together = results_folder + os.path.sep + 'together' + os.path.sep + f'{endpoint}'
                create_dir_for_io(merged_files_folder_selcols_together)

                merged_output_selcols.to_csv(merged_files_folder_selcols_together + os.path.sep + file_name_output + f'_{complete_endpoint}.csv', sep = ';', index = False)





                if merge_ONTOX == True:
                    processed_dfs[tag_tool] = {}

                    processed_dfs[tag_tool]['df'] = merged_output_selcols
                    processed_dfs[tag_tool]['pred_units'] = predicted_units
                    processed_dfs[tag_tool]['des_units'] = desired_units
                    
                if merge_SITOLUB == True:
                    processed_dfs[tag_tool] = {}

                    processed_dfs[tag_tool]['df'] = merged_output_selcols
                    processed_dfs[tag_tool]['pred_units'] = predicted_units
                    processed_dfs[tag_tool]['des_units'] = desired_units    
                    
                # if f'ID_{tag_tool}' in merged_output_selcols.columns:
                #     processed_dfs[tag_tool]['idindf'] = 'yes'
                # else:
                #     processed_dfs[tag_tool]['idindf'] = 'no'
        
        elif tool=="VEGA_tool":
            tag_tool = tool.split('_')[0]
            equivalence_in_tool = df_tools_enpoint_names.loc[tag,tool]
            tool_folder = results_folder + os.path.sep + tool
            files_to_predict_folder = tool_folder + os.path.sep + 'inputs'
            files_predicted_folder = tool_folder + os.path.sep + 'predictions'
            #collect files to be predicted
            files_topredict_property = collect_files(files_to_predict_folder)



            #collect files predicted
            files_predicted_property = collect_files(files_predicted_folder)
            # print('patata')
            merged_input, file_name_input, models = merge_files(files_topredict_property, files_to_predict_folder, tool, equivalence_in_tool, mode = 'inputfiles')
                # merge output files

            merged_output, file_name_output, models = merge_files(files_predicted_property, files_predicted_folder, tool, equivalence_in_tool, mode = 'outputfiles')
            
            
            if equivalence_in_tool=="MM":
                # cols=merged_output.columns
                # print("AAAAAAAAAAAAAAAAAAAAAAA",merged_output.columns)

                merged_output.rename(columns={"SMILES":f"SMILES_{tag_tool}","tId":f"ID_{tag_tool}"},inplace=True)
                merged_output["SMILES"]=merged_output[f"SMILES_{tag_tool}"]
                merged_output["SMILES"]=merged_output["SMILES"].apply(sanitize_smiles)
                # print("AAAAAAAAAAAAAAAAAAAAAAA",merged_output.columns)

                for model2 in models:
                    col_exp = f'Exp_{model2}-assessment'
                    tagnuevo_exp = f'exp_{model2}_{tag_tool}_{endpoint}'

                    col_AD = f'AD_{model2}-assessment'
                    tagnuevo_AD = f'AD_{model2}_{tag_tool}_{endpoint}'

                    col_pred = f'Pred_{model2}-prediction'
                    tagnuevo_pred = f'pred_{model2}_{tag_tool}_{endpoint}'

                    merged_output.rename(columns={
                        col_exp: tagnuevo_exp,
                        col_AD: tagnuevo_AD,
                        col_pred: tagnuevo_pred
                    }, inplace=True)
                # print(merged_output.columns)
                merged_files_folder = results_folder + os.path.sep + tool + os.path.sep + 'merged_predictions_all'
                create_dir_for_io(merged_files_folder)
                merged_output.to_csv(merged_files_folder + os.path.sep + file_name_output + '_allcolumns.csv', sep = ';', index = False)


                if merge_SITOLUB == True:
                    processed_dfs[tag_tool] = {}

                    processed_dfs[tag_tool]['df'] = merged_output
                    processed_dfs[tag_tool]['pred_units'] = np.nan
                    processed_dfs[tag_tool]['des_units'] = np.nan  

                
        elif tool == "ECOSAR_tool":
            print("a")
        
        elif tool=="EXPERIMENTAL_tool":
            # print("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
            tag_tool = tool.split('_')[0]
            
            equivalence_in_tool = df_tools_enpoint_names.loc[tag,tool]
            
            
            tool_folder = results_folder + os.path.sep + tool
            files_to_predict_folder = tool_folder + os.path.sep + 'inputs'
            files_predicted_folder = tool_folder + os.path.sep + 'predictions'
            #collect files to be predicted
            files_topredict_property = collect_files(files_to_predict_folder)



            #collect files predicted
            files_predicted_property = collect_files(files_predicted_folder)
            # print('patata')
            # merged_input, file_name_input, models = merge_files(files_topredict_property, files_to_predict_folder, tool, equivalence_in_tool, mode = 'inputfiles')
                # merge output files

            
            
            
            merged_output, file_name_output, models = merge_files(files_predicted_property, files_predicted_folder, tool, equivalence_in_tool, mode = 'outputfiles')
            # print("AAAAAAAAAAAAAAAAAAAAAAA",merged_output.columns)
            if equivalence_in_tool=="EXPERIMENTAL":
                # print(merged_output.columns)
                merged_output.rename(columns={"SMILES":f"SMILES_{tag_tool}","ID":f"ID_{tag_tool}"},inplace=True)
                merged_output["SMILES"]=merged_output[f"SMILES_{tag_tool}"].apply(sanitize_smiles)
                # print("AAAAAAAAAAAAAAAAAAAAAAA",merged_output)
                merged_files_folder = results_folder + os.path.sep + tool + os.path.sep + 'merged_predictions_all'
                create_dir_for_io(merged_files_folder)
                merged_output.to_csv(merged_files_folder + os.path.sep + file_name_output + '_allcolumns.csv', sep = ';', index = False)


                if merge_SITOLUB == True:
                    processed_dfs[tag_tool] = {}

                    processed_dfs[tag_tool]['df'] = merged_output
                    processed_dfs[tag_tool]['pred_units'] = np.nan
                    processed_dfs[tag_tool]['des_units'] = np.nan  

    
    #     if merge_SITOLUB == True:
    #         validate_df[tag_tool] = {}

    #         validate_df[tag_tool]['df'] = merged_output

    #         # print('Pruebas JL')
    #         # print(tag_tool)
    #         # print(validate_df.keys())
            
    #         # processed_dfs[tag_tool]['pred_units'] = predicted_units
    #         # processed_dfs[tag_tool]['des_units'] = desired_units
    #     # if f'ID_{tag_tool}' in merged_output_selcols.columns:
    #     #     processed_dfs[tag_tool]['idindf'] = 'yes'
    #     # else:
    #     #     processed_dfs[tag_tool]['idindf'] = 'no'
    
    #     if not val.empty:
    #         # print(val)
    #         merged_output=merged_output.drop(columns=['index'], errors='ignore')
           
            
    #         val = pd.concat([val, merged_output], axis=1)
    #         val.to_csv("a.csv",sep=";",encoding="utf8",index=False)
    #     else:
    #         val = merged_output
    #         val.set_index('SMILES')
        
    # # print(val)
    # # print(validate_df.keys())

    # print('\tAll files processed!!!!!')


#%% manage_for_ONTOX
    # if merge_SITOLUB == True:
        
    #     print('ahorasiquemehallo')
    #     print(' \n\n [++] Creating validation file (SITOLUB)')
        
    #     matched_files_folder = results_folder + os.path.sep + 'together' + os.path.sep + 'matched'
    #     create_dir_for_io(matched_files_folder)
    #     print(f'\t\t{matched_files_folder}')
    #     original_file_name = config_df.loc[dataset, 'file']
    #     print(f'\t\t{original_file_name}')


    #     original_file = pd.read_csv(initial_data_folder + os.path.sep + original_file_name, sep = ';', low_memory=False)
    #     print(original_file)
    #     print(os.scandir(initial_data_folder),"FFFFFFFFFFFFFFFFF")


    #     df_duplis = original_file.loc[original_file['SMILES'].duplicated(keep=False)]
    #     if df_duplis.shape[0] != 0:
    #         print('\t\t WARNING!!!! you have duplicated molecules in your prediction files ')
    #         print('\t\t processing cannot continue')

    #     else:
    #         print('estoy bien ')


    #     all_dfs_with_ids  = [original_file]
    #     # print(all_dfs_with_ids)
    #     all_dfs_wo_ids  = []
    #     print(all_dfs_wo_ids)

    #     for key, value in processed_dfs.items():
    #         print(key,value)
    #         print(f'\t\t[++] Processing {key} for ONTOX')

    #         df = value['df']
    #         # print(df)



    #         final_cols = [f'SMILES_{key}', f'pred_{key}_{endpoint}_converted']

    #         if f'AD_{key}_{endpoint}' in list(df.columns):


    #             final_cols.append(f'AD_{key}_{endpoint}')

    #         if f'in_TS_{key}_{endpoint}' in list(df.columns):
    #             final_cols.append(f'in_TS_{key}_{endpoint}')

    #         if f'ID_{key}' in list(df.columns):


    #             final_cols = [f'ID_{key}'] + final_cols
    #             # print(final_cols)
    #             # print(list(df.columns))
    #             # print(final_cols)


    #             df_to_merge = df[final_cols]

    #             df_to_merge.rename(columns = {
    #                         f'pred_{key}_{endpoint}_converted': f'{key}_{endpoint}_pred',
    #                         f'AD_{key}_{endpoint}': f'{key}_{endpoint}_AD',
    #                         f'in_TS_{key}_{endpoint}': f'{key}_{endpoint}_TS',

    #                         }, inplace = True)

    #             df_to_merge = df_to_merge.set_index(f'ID_{key}')

    #             df_duplis = df_to_merge.loc[df_to_merge.index.duplicated(keep=False)]
    #             if df_duplis.shape[0] != 0:
    #                 print('\t\t WARNING!!!! you have duplicated molecules in your prediction files ')
    #                 print('\t\t processing cannot continue')

    #             # else:
    #             #     print('estoy bien ')

    #             all_dfs_with_ids.append(df_to_merge)




    #             original_file.set_index('ID_all', inplace = True)
    #             original_file = pd.concat([original_file,df_to_merge], axis = 1, join="outer")
    #             original_file.index.name = 'ID_all'
    #             original_file.reset_index(inplace = True)

    #         else:

    #             # print(list(df.columns))
    #             # print(final_cols)




    #             df_to_merge = df[final_cols]


    #             df_to_merge.rename(columns = {
    #                                 f'pred_{key}_{endpoint}_converted': f'{key}_{endpoint}_pred',
    #                                 f'AD_{key}_{endpoint}': f'{key}_{endpoint}_AD',
    #                                 f'in_TS_{key}_{endpoint}': f'{key}_{endpoint}_TS',

    #                     }, inplace = True)

    #             df_to_merge['SMILES'] = df_to_merge[f'SMILES_{key}']
    #             # df_to_merge['SMILES'] = df_to_merge.apply(lambda row: sanitize_smiles(row[f'SMILES_{key}']), axis=1)



    #             # df_to_merge.index.name = f'ID_{key}'

    #             # df_to_merge.reset_index(drop = False, inplace = True)

    #             df_to_merge.set_index('SMILES', inplace = True)

    #             # manage duplicated molecules



    #             df_duplis_onerepe = df_to_merge.loc[df_to_merge.index.duplicated(keep='first')]
    #             df_duplis = df_to_merge.loc[df_to_merge.index.duplicated(keep=False)]


    #             if df_duplis.shape[0] != 0:
    #                 print('\t\t WARNING!!!! you have duplicated molecules in your prediction files ')
    #                 print(f'\t\t a total of {df_duplis_onerepe.shape[0]} have been skeeped. A new file has been created')

    #                 df_duplis.to_csv(matched_files_folder + os.path.sep + f'{tag}_{key}_skeepedbyduplication.csv', sep = ';', index = False)


    #                 df_to_merge = df_to_merge.loc[~df_to_merge.index.duplicated(keep='first')]


    #             # manage molecules where smiles are changes by admetlab

    #             diff = df_to_merge[~df_to_merge[f'SMILES_{key}'].isin(original_file['SMILES'])]

    #             if diff.shape[0] != 0:
    #                 print('\t\t WARNING!!!! some molecules cannot be matched ')
    #                 print(f'\t\t a total of {diff.shape[0]} have been skeeped. A new file has been created')
    #                 diff.to_csv(matched_files_folder + os.path.sep + f'{tag}_{key}_skeepednotmatched.csv', sep = ';', index = True)



    #             onlyequal = df_to_merge[df_to_merge[f'SMILES_{key}'].isin(original_file['SMILES'])]

    #             original_file.set_index('SMILES', inplace = True)

    #             all_dfs_wo_ids.append(onlyequal)


    #             # onlyequal.drop(f'ID_{key}', axis = 1)


    #             original_file = pd.concat([original_file,onlyequal], axis = 1, join="outer")

    #             original_file.index.name = 'SMILES'
    #             original_file.reset_index(inplace = True)
    #             original_file.set_index('ID_all', inplace = True)
    #             original_file.reset_index(inplace = True)






    #         original_file.to_csv(matched_files_folder + os.path.sep + f'{tag}.csv', sep = ';', index = False)
    
    
    if merge_ONTOX == True or merge_SITOLUB ==True:

        print(' \n\n [++] Matching with original file (ONTOX)')

        matched_files_folder = results_folder + os.path.sep + 'together' + os.path.sep + 'matched'
        create_dir_for_io(matched_files_folder)

        original_file_name = config_df.loc[dataset, 'file']
        print(f'\t\t{original_file_name}')


        original_file = pd.read_csv(initial_data_folder + os.path.sep + original_file_name, sep = ';', low_memory=False)



        df_duplis = original_file.loc[original_file['SMILES'].duplicated(keep=False)]
        if df_duplis.shape[0] != 0:
            print('\t\t WARNING!!!! you have duplicated molecules in your prediction files ')
            print('\t\t processing cannot continue')

        # else:
        #     print('estoy bien ')


        all_dfs_with_ids  = [original_file]

        all_dfs_wo_ids  = []
        
        
         
        


        for key, value in processed_dfs.items():

            print(f'\t\t[++] Processing {key} for ONTOX')

            df = value['df']
            
            
            if key == 'VEGA' :
                
                print('arbol')
                
                # original_file = pd.concat([original_file,onlyequal], axis = 1, join="outer")
    
                # original_file.index.name = 'SMILES'
                # original_file.reset_index(inplace = True)
                # original_file.set_index('ID_all', inplace = True)
                # original_file.reset_index(inplace = True)
                
                
                
            else:    


    
    
                final_cols = [f'SMILES_{key}', f'pred_{key}_{endpoint}_converted']
    
                if f'AD_{key}_{endpoint}' in list(df.columns):
    
    
                    final_cols.append(f'AD_{key}_{endpoint}')
    
                if f'in_TS_{key}_{endpoint}' in list(df.columns):
                    final_cols.append(f'in_TS_{key}_{endpoint}')
    
                if f'ID_{key}' in list(df.columns):
    
    
                    final_cols = [f'ID_{key}'] + final_cols
                    # print(list(df.columns))
                    # print(final_cols)
    
    
                    df_to_merge = df[final_cols]
    
                    df_to_merge.rename(columns = {
                                f'pred_{key}_{endpoint}_converted': f'{key}_{endpoint}_pred',
                                f'AD_{key}_{endpoint}': f'{key}_{endpoint}_AD',
                                f'in_TS_{key}_{endpoint}': f'{key}_{endpoint}_TS',
    
                                }, inplace = True)
    
                    df_to_merge = df_to_merge.set_index(f'ID_{key}')
    
                    df_duplis = df_to_merge.loc[df_to_merge.index.duplicated(keep=False)]
                    if df_duplis.shape[0] != 0:
                        print('\t\t WARNING!!!! you have duplicated molecules in your prediction files ')
                        print('\t\t processing cannot continue')
    
                    # else:
                    #     print('estoy bien ')
    
                    all_dfs_with_ids.append(df_to_merge)
    
    
    
    
                    original_file.set_index('ID_all', inplace = True)
                    original_file = pd.concat([original_file,df_to_merge], axis = 1, join="outer")
                    original_file.index.name = 'ID_all'
                    original_file.reset_index(inplace = True)
    
                else:
    
                    # print(list(df.columns))
                    # print(final_cols)
    
    
    
    
                    df_to_merge = df[final_cols]
    
    
                    df_to_merge.rename(columns = {
                                        f'pred_{key}_{endpoint}_converted': f'{key}_{endpoint}_pred',
                                        f'AD_{key}_{endpoint}': f'{key}_{endpoint}_AD',
                                        f'in_TS_{key}_{endpoint}': f'{key}_{endpoint}_TS',
    
                            }, inplace = True)
    
                    df_to_merge['SMILES'] = df_to_merge[f'SMILES_{key}']
                    # df_to_merge['SMILES'] = df_to_merge.apply(lambda row: sanitize_smiles(row[f'SMILES_{key}']), axis=1)
    
    
    
                    # df_to_merge.index.name = f'ID_{key}'
    
                    # df_to_merge.reset_index(drop = False, inplace = True)
    
                    df_to_merge.set_index('SMILES', inplace = True)
    
                    # manage duplicated molecules
    
    
    
                    df_duplis_onerepe = df_to_merge.loc[df_to_merge.index.duplicated(keep='first')]
                    df_duplis = df_to_merge.loc[df_to_merge.index.duplicated(keep=False)]
    
    
                    if df_duplis.shape[0] != 0:
                        print('\t\t WARNING!!!! you have duplicated molecules in your prediction files ')
                        print(f'\t\t a total of {df_duplis_onerepe.shape[0]} have been skeeped. A new file has been created')
    
                        df_duplis.to_csv(matched_files_folder + os.path.sep + f'{tag}_{key}_skeepedbyduplication.csv', sep = ';', index = False)
    
    
                        df_to_merge = df_to_merge.loc[~df_to_merge.index.duplicated(keep='first')]
    
    
                    # manage molecules where smiles are changes by admetlab
    
                    diff = df_to_merge[~df_to_merge[f'SMILES_{key}'].isin(original_file['SMILES'])]
    
                    if diff.shape[0] != 0:
                        print('\t\t WARNING!!!! some molecules cannot be matched ')
                        print(f'\t\t a total of {diff.shape[0]} have been skeeped. A new file has been created')
                        diff.to_csv(matched_files_folder + os.path.sep + f'{tag}_{key}_skeepednotmatched.csv', sep = ';', index = True)
    
    
    
                    onlyequal = df_to_merge[df_to_merge[f'SMILES_{key}'].isin(original_file['SMILES'])]
    
                    original_file.set_index('SMILES', inplace = True)
    
                    all_dfs_wo_ids.append(onlyequal)
    
    
                    # onlyequal.drop(f'ID_{key}', axis = 1)
    
    
                    original_file = pd.concat([original_file,onlyequal], axis = 1, join="outer")
    
                    original_file.index.name = 'SMILES'
                    original_file.reset_index(inplace = True)
                    original_file.set_index('ID_all', inplace = True)
                    original_file.reset_index(inplace = True)
    
    




            original_file.to_csv(matched_files_folder + os.path.sep + f'{tag}.csv', sep = ';', index = False)

#%%


        # for dftomerge in all_dfs_with_ids:



        #     df_duplis = dftomerge.loc[df.index.duplicated(keep=False)]

        #     if df_duplis.shape[0] != 0:
        #         print('WARNING!!!! you have duplicated molecules in your prediction files ')


#%%



        # all_merged = pd.concat(all_dfs_with_ids, axis = 1, join="outer")

        # all_merged.index.name = 'ID_all'

        # all_merged.reset_index(inplace = True)

        # all_merged_bysmiles = all_merged.set_index('SMILES')

        # to_final_merge = [all_merged_bysmiles] + all_dfs_wo_ids

        # all_merged_bysmiles_merged = pd.concat(to_final_merge, axis = 1, join="outer")        # print(all_merged.shape)

        # matched_files_folder = results_folder + os.path.sep + 'together' + os.path.sep + 'matched'
        # all_merged_bysmiles_merged.to_csv(matched_files_folder + os.path.sep + f'{tag}.csv', sep = ';', index = False)






#%%
















