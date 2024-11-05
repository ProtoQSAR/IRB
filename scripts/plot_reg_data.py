# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 08:06:19 2024

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
import matplotlib.pyplot as plt
''' requires to install pip install xlrd '''
##############################################################################
################################ INITIAL VARIABLES ###########################
parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

input_path = './compare_ONGOING/results'

df = pd.read_csv(input_path + os.path.sep + 'all_merged_TOX_MRDD.csv', sep = ';')

df_valid = df[df['valid/outlier'] == 'valid']

plt.hist(df_valid['y_consensus'])

df_valid_sorted = df_valid.sort_values('y_consensus',ascending=True)

df_valid_sorted.reset_index(inplace = True, drop = True)

plt.scatter(df_valid_sorted.index, df_valid_sorted['y_consensus']) 

plt.show()



df_valid_deepk = df_valid[['ID_all','y_DeepPK_nointraoutliers']]

df_valid_deepk_notnull = df_valid_deepk[pd.notnull(df_valid_deepk['y_DeepPK_nointraoutliers'])]

plt.hist(df_valid_deepk_notnull['y_DeepPK_nointraoutliers'])

df_valid_deepk_notnull_sorted = df_valid_deepk_notnull.sort_values('y_DeepPK_nointraoutliers',ascending=True)

df_valid_deepk_notnull_sorted.reset_index(inplace = True, drop = True)

plt.scatter(df_valid_deepk_notnull_sorted.index, df_valid_deepk_notnull_sorted['y_DeepPK_nointraoutliers']) 

plt.show() 




df_valid_DSSTOX = df_valid[['ID_all','y_DSSTOX_nointraoutliers']]

df_valid_DSSTOX_notnull = df_valid_DSSTOX[pd.notnull(df_valid_DSSTOX['y_DSSTOX_nointraoutliers'])]

plt.hist(df_valid_DSSTOX_notnull['y_DSSTOX_nointraoutliers'])

df_valid_DSSTOX_notnull_sorted = df_valid_DSSTOX_notnull.sort_values('y_DSSTOX_nointraoutliers',ascending=True)

df_valid_DSSTOX_notnull_sorted.reset_index(inplace = True, drop = True)

plt.scatter(df_valid_DSSTOX_notnull_sorted.index, df_valid_DSSTOX_notnull_sorted['y_DSSTOX_nointraoutliers']) 

plt.show() 