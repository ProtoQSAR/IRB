# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:13:33 2024

This script will merge files for train, test and validation from DeepPK database

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

import numpy as np
import pandas as pd
import time
''' requires to install pip install xlrd '''
##############################################################################
################################ INITIAL VARIABLES ###########################
parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')

inputpath = '..' + os.path.sep + 'datasets_for_modelling'
outputpath =  'pre_preprocessed_ONGOING' + os.path.sep + 'input_data'

import pandas as pd


#properties = ['fdamdd_reg', 'cyp1a2_inhibitor','cyp2c9_inhibitor', 'cyp2c19_inhibitor', 'cyp2d6_inhibitor', 'cyp3a4_inhibitor', 'herg']

properties = ['nr_ar'] #for HYPIEND

for properti in properties:
    df_train = pd.read_csv(inputpath + os.path.sep + f'{properti}_train.csv', sep = ',')
    df_test = pd.read_csv(inputpath + os.path.sep + f'{properti}_test.csv', sep = ',')
    df_val = pd.read_csv(inputpath + os.path.sep + f'{properti}_val.csv', sep = ',')
    
    df_all = pd.concat([df_train,df_test,df_val], axis = 0)
    
    
    if properti != 'fdamdd_reg': #because it is regression data
        print(properti)
        print('\t\t negative', df_all['label'].value_counts()[0])
        print('\t\t positive', df_all['label'].value_counts()[1])
    
    df_all.index.name = 'ID'
    
    df_all.to_csv(outputpath + os.path.sep + f'DeepPK_{properti}_all.csv', sep = ';', index= True)