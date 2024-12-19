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
from io import StringIO
import re
from ast import literal_eval
''' requires to install pip install xlrd '''
##############################################################################

################################ INITIAL VARIABLES ###########################
parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

data_folder =  '.' + os.path.sep + 'files'

mixtures_folder =  '.' + os.path.sep + 'pre_processing'
cluster_folder =  '..' + os.path.sep + 'clustering' + os.path.sep + 'clusters_kmeans'

predictions_folder =  '..' + os.path.sep + 'scripts' + os.path.sep + 'reimputation_and_prediction' + os.path.sep + 'predictions' + os.path.sep + 'df_predicted'

results_folder_with_clustering = '.' + os.path.sep + 'clusterized'

name = 'IRB_library_'
###############################  FUNCTIONS ###################################

def retrieve_smiles(mol):

    return Chem.MolToSmiles(mol)


def create_folder(folder):

    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Output directory created: {folder}")
    else:
        print(f"Output directory already exists: {folder}")


def getidentifier(smiles):
    
    identif = df_clustering_asdict_smiles[smiles]['ID NUMBER']
    
    return identif

##############################################################################


create_folder(results_folder_with_clustering)


#%%
##############################################################################
############################ READ INITIAL SDF FILES ##########################
##############################################################################

print('\nREADING INITIAL SDF FILES')

############################## file 1 processing #############################

file1 = '2021 LIB_47489 CMPDS_26092024.sdf'

print(f'\n[+] Opening "{file1}" file')

df1 = PandasTools.LoadSDF(data_folder + os.path.sep + file1)

df1['SMILES']  = df1.apply(lambda row: retrieve_smiles(row['ROMol']), axis=1)

df1.insert(0, 'ID_SET', ['set_1']*df1.shape[0])

df1.index.name = 'ID_SET_COMPOUND'

output_file1_name = file1.split('.')[0]
df1.to_csv(results_folder_with_clustering+ os.path.sep + name + output_file1_name + '_tobeclustered.csv', sep = ';')

print(f'\t[++] {file1}: {df1.shape[0]} compounds')

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
df2.to_csv(results_folder_with_clustering+ os.path.sep + name + output_file2_name + '_tobeclustered.csv', sep = ';')

print(f'\t[++] {file2}: {df2.shape[0]} compounds')

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
df3.to_csv(results_folder_with_clustering+ os.path.sep + name + output_file3_name + '_tobeclustered.csv', sep = ';')

print(f'\t[++] {file3}: {df3.shape[0]} compounds')
#%%
##############################################################################

############################## file 4 processing #############################

# this file has some incorrect molecules, thus, it uses different approach to be load in the final step

file4 = 'focus 10 mM_10665 cmpds_29092024_sent.sdf'

print(f'\n[+] Opening "{file4}" file')

df4 = PandasTools.LoadSDF(data_folder + os.path.sep + file4)

df4.rename(columns = {'IRB plate ': 'IRB PLATE', 'Library': 'LIBRARY'}, inplace = True)

df4['SMILES']  = df4.apply(lambda row: retrieve_smiles(row['ROMol']), axis=1)

df4.insert(0, 'ID_SET', ['set_4']*df4.shape[0])

df4.index.name = 'ID_SET_COMPOUND'


output_file4_name = file4.split('.')[0]
df4.to_csv(results_folder_with_clustering+ os.path.sep + name + output_file4_name + '_tobeclustered.csv', sep = ';')

print(f'\t[++] {file4}: {df4.shape[0]} compounds')

#%%
df4_withincorrect = PandasTools.LoadSDF(data_folder + os.path.sep + file4, molColName = None)

df4_incorrect = df4_withincorrect[~df4_withincorrect['ID NUMBER'].isin(df4['ID NUMBER'])]



df4_incorrect.to_csv(results_folder_with_clustering+ os.path.sep + name + output_file4_name + '_molerror.csv', sep = ';')



##############################################################################


#%%
##############################################################################
######################### INCORPORATE CLUSTERING INFO ########################
##############################################################################
print('\nINCORPORATING CLUSTERING INFO')

print('[+] Retrieving clustering info')

# df_clustering = pd.read_excel(cluster_folder +  os.path.sep + 'IRB_all_clustering40K_sorted.xlsx')
df_clustering = pd.read_csv(cluster_folder +  os.path.sep + 'IRB_library_merged_clustering40K_sorted.csv', sep=';')

# df_clustering_short = df_clustering.iloc[1660:1670,:]

dfclustering_dict = {}

for idx, row in df_clustering.iterrows():
    
    try:
        identifier = row["ID NUMBER"]
        splitted = [i for i in literal_eval(identifier)]
        for ids in splitted:
            dfclustering_dict[ids] = {}
            dfclustering_dict[ids]['SMILES_new'] = row['SMILES']
            dfclustering_dict[ids]['Cluster'] = row['Cluster']
            dfclustering_dict[ids]['is_centroid'] = row['is_centroid']
        
    except:
        ids = row["ID NUMBER"]
        dfclustering_dict[ids] = {}
        dfclustering_dict[ids]['SMILES_new'] = row['SMILES']
        dfclustering_dict[ids]['Cluster'] = row['Cluster']
        dfclustering_dict[ids]['is_centroid'] = row['is_centroid']




mixtures_identified = pd.read_csv(mixtures_folder + os.path.sep + 'IRB_library_merged_all-mixtures_for_manual_analysis.csv', sep = ';')

mixtures_identified_dic = mixtures_identified.set_index('ID NUMBER').T.to_dict(orient = 'series')

df1_dict = df1.T.to_dict()
df2_dict = df2.T.to_dict()
df3_dict = df3.T.to_dict()
df4_dict = df4.T.to_dict()


print('[+] Including clustering info in original files')

dictios_names = ['2021 LIB_47489 CMPDS_26092024_clustered', 'all new library_104017 CMPDS_26092024_clustered', 'focus 2 mM_70 cmpds_29092024_sent_clustered', 'focus 10 mM_10665 cmpds_29092024_sent_clustered']

dictios = [df1_dict,df2_dict,df3_dict,df4_dict]

dictios_with_cluster = []
# dictios_names = ['focus 2 mM_70 cmpds_29092024_sent_clustered']

# dictios = [df3_dict]



for name, dictio in zip(dictios_names, dictios):
    
    print(f'\t[++] Working on {name} sdf')

    for idx2, row2 in dictio.items():
        if row2['ID NUMBER'] in dfclustering_dict.keys():
            flag = '-'
            smiles_new = dfclustering_dict[row2['ID NUMBER']]['SMILES_new']
            cluster = dfclustering_dict[row2['ID NUMBER']]['Cluster']
            centroid = dfclustering_dict[row2['ID NUMBER']]['is_centroid']
            
            
        else:
            if row2['ID NUMBER'] in mixtures_identified_dic.keys():
                flag = 'mixture'
                smiles_new = '-'
                cluster = '-'
                centroid = '-'
            else:
                if row2['ID'] =='empty mol':
                    flag = 'empty mol'
                    smiles_new = '-'
                    cluster = '-'
                    centroid = '-'
                else:
                    flag = 'inorganic_organometallic_invalid'
                    smiles_new = '-'
                    cluster = '-'
                    centroid = '-'
                    
        dictio[idx2]['SMILES_new'] = smiles_new
        dictio[idx2]['Cluster'] = cluster
        dictio[idx2]['is_centroid'] = centroid
        dictio[idx2]['flag'] = flag
        
        dictio[idx2]['Clustering_method'] = 'Clustering performed by ProtoQSAR, S.L., 2024'
        
    dictios_with_cluster.append(dictio)

    df_again = pd.DataFrame.from_dict(dictio).T
    
    print(f'\t\t {df_again.shape[0]} compounds')
    
    
    strip_number = re.compile(r"^(>  <[\w\s]+>)(\s+\(\d+\)\s*)$", re.MULTILINE)
    
    with StringIO() as buf:
        PandasTools.WriteSDF(df_again, buf, properties=list(df_again.columns))
        sdf = buf.getvalue()
    sdf = strip_number.sub(r"\1", sdf)
    with open(results_folder_with_clustering + os.path.sep + f'{name}_clustered.sdf', "w") as hnd:
        hnd.write(sdf)
    

#%%
##############################################################################
########################### INCORPORATE PREDICTIONS ##########################
##############################################################################
print('\nINCORPORATING PREDICTIONS')

print('[+] Retrieving, merging predictions and assigning identifier')

df_clustering_asdict_smiles = df_clustering.set_index('SMILES').T.to_dict(orient = 'dict') # to retrieve identifiers


for root_pred, dirs_pred, files_pred in os.walk(predictions_folder):
    continue

predictions_parts = [x.split('-')[0] for x in files_pred]

set_predictions_parts = set(predictions_parts)

endpoints = [x.split('-')[1] for x in files_pred]

set_endpoints_aslist = list(set(endpoints))

#%%

info_dict = {'TK_F20': {'explanation': 'Bioavailability 20% [binary: (1:positive, 0: negative), positive if Bioavailability > 20%]'},
             'TK_F30': {'explanation': 'Bioavailability 30% [binary: (1:positive, 0: negative),  positive if Bioavailability > 30%]'},
             'TK_FU': {'explanation': 'Fraction unbound (FU) to plasma proteins [ratio, dimensionless]'},

             'TK_HLM': {'explanation': 'Human Liver Microsomal Stability [binary: (1:stable, 0: unstable)], stable if CLint < 20 mL/min/kg]'},
             'TK_logKp': {'explanation': 'Skin permeability [log(cm/h)]'},
             'TK_Caco2': {'explanation': ' Caco-2 permeability [log(cm/s)]'},
             'TK_VDss': {'explanation': 'Volumne of distribution [log(L/kg)]'},
             'TK_HIA': {'explanation': 'Human intestinal absorption [binary: (1:positive, 0: negative), positive if HIA% > 30%]'},
             'TK_BBB': {'explanation': ' Blood-brain barrier penetration [binary: (1:positive, 0: negative), positive if LogBBB ≥ -1]'},
             
             'TK_CYP2C19inh': {'explanation': 'Cytochrome P450 2C19 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},
             'TK_CYP1A2inh': {'explanation': 'Cytochrome P450 1A2 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},             
             'TK_CYP2D6inh': {'explanation': 'Cytochrome P450 2D6 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},            
             'TK_CYP3A4inh': {'explanation': 'Cytochrome P450 3A4 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},          
             'TK_CYP2C9inh': {'explanation': 'Cytochrome P450 2C9 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},

             'TK_CYP2C19sub': {'explanation': 'Cytochrome P450 2C19 substrate [binary: (1:substrate, 0: non-substrate)]'},
             'TK_CYP1A2sub': {'explanation': 'Cytochrome P450 1A2 substrate [binary: (1:substrate, 0: non-substrate)]'},
             'TK_CYP2D6sub': {'explanation': 'Cytochrome P450 2D6 substrate [binary: (1:substrate, 0: non-substrate)]'},
             'TK_CYP3A4sub': {'explanation': 'Cytochrome P450 3A4 substrate [binary: (1:substrate, 0: non-substrate)]'},
             'TK_CYP2C9sub': {'explanation': 'Cytochrome P450 2C9 substrate [binary: (1:substrate, 0: non-substrate)]'},

             'TK_OATP1B1inh': {'explanation': 'Organic anion transporting polypeptides 1B1 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},
             'TK_OATP1B3inh': {'explanation': 'Organic anion transporting polypeptides 1B3 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},

             'TK_Pgpsub': {'explanation': 'P-glycoprotein substrate [binary: (1:substrate, 0: non-substrate)]'}, 
             'TK_Pgpinh': {'explanation': 'P-glycoprotein inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},
             'TK_BSEPinh': {'explanation': 'Bile salt export pump substrate [binary: (1:substrate, 0: non-substrate)]'},
             

             'TOX_MRDD': {'explanation': 'Maximum Recommended Daily Dose [log(mg/kg/day)]'},
             'TOX_hERGinh': {'explanation': 'Human ether-a-go-go-related gene (hERG) inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},
             'TOX_Cav12inh': {'explanation': 'voltage-gated calcium channel Ca V 1.2 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'},             
             'TOX_Nav15inh': {'explanation': 'voltage-gated sodium channel Na(v) 1.5 inhibitor [binary: (1:inhibitor, 0: non-inhibitor)]'}
             
             }


#%%


for endpoint in set_endpoints_aslist:

    
    if endpoint == 'VDss_log10_sin_outliers':
        clean_endpoint = 'TK_VDss'
    else:
        clean_endpoint_list = endpoint.split('_')
        clean_endpoint = clean_endpoint_list [0] + '_' + clean_endpoint_list [1]
    

    
    info_dict[clean_endpoint]['original_name'] = endpoint



#%%


#%%
merged_by_endpoints = []

mergedallendpoints = pd.DataFrame()
    
for cleanendpoint, value in info_dict.items():
    
    print(f'\t[++] Merging predictions for {endpoint}')
    
    endpoint = info_dict[cleanendpoint]['original_name']
    
    mergedbyendpoint = pd.DataFrame()
    
    for part in set_predictions_parts:

        
        df_pred = pd.read_csv(predictions_folder + os.path.sep + f'{part}-{endpoint}-predictedAD.csv', sep = ';')
        
        sel_cols = ['SMILES', f'{endpoint} experimental value', f'{endpoint} predicted value', 'AD_ProtoPRED']
        
        df_pred_selcols = df_pred[sel_cols]
        
        df_pred_selcols.rename(columns = {f'{endpoint} experimental value': f'{cleanendpoint}_experimental_value',
                                          f'{endpoint} predicted value' : f'{cleanendpoint}_predicted_value',
                                          'AD_ProtoPRED': f'{cleanendpoint}_AD'}, inplace = True)
        
        #df_pred_selcols[f'{cleanendpoint}_info'] = info_dict[cleanendpoint]['explanation'] * df_pred_selcols.shape[0] 

        mergedbyendpoint = pd.concat([mergedbyendpoint,df_pred_selcols], axis = 0)
        
    print(f'\t\t {mergedbyendpoint.shape[0]} compounds')
    
    mergedbyendpoint.set_index('SMILES', inplace = True)
    
    mergedallendpoints = pd.concat([mergedallendpoints,mergedbyendpoint], axis = 1)
    
    print(mergedallendpoints.shape)
    
    
dict_allpredictions = mergedallendpoints.T.to_dict(orient = 'dict') 
    
#%%  


print('[+] Include prediction on sdf')

for name, dictiocluster in zip(dictios_names, dictios_with_cluster):
    
    print(f'\t[++] Working on {name} sdf')


    for idx3, row3 in dictiocluster.items():
        if row3['SMILES_new'] in dict_allpredictions.keys():

            for key1, value1 in dict_allpredictions[row3['SMILES_new']].items():
                dictiocluster[idx3][key1] = value1
            

        else:
            
            if row3['flag'] != '-':
                continue
            else:
                print('SOS')
                
        dictiocluster[idx3]['ADME_prediction_method'] = 'ADMET prediction performed with ProtoPRED (R) v1.0'
                
                
    df_again2 = pd.DataFrame.from_dict(dictiocluster).T


            
    print(f'\t\t {df_again2.shape[0]} compounds')
    
    
    strip_number = re.compile(r"^(>  <[\w\s]+>)(\s+\(\d+\)\s*)$", re.MULTILINE)
    
    with StringIO() as buf:
        PandasTools.WriteSDF(df_again2, buf, properties=list(df_again2.columns))
        sdf = buf.getvalue()
    sdf = strip_number.sub(r"\1", sdf)
    with open(results_folder_with_clustering + os.path.sep + f'{name}_predicted_renamed.sdf', "w") as hnd:
        hnd.write(sdf)            


##############END OF THE SCRIPT###############################################







